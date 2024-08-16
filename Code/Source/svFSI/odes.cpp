#include "odes.h"

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if (from.empty())
        return;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

Expression create_symbolic_expression(const std::string& expr_str, const std::unordered_map<std::string, RCP<const Basic>>& symbols_map) {
    std::string parsed_expr = expr_str;
    for (const auto& pair : symbols_map) {
        std::string symbol_name = pair.first;
        std::string symbol_rep = SymEngine::str(*pair.second);
        size_t pos = 0;
        while ((pos = parsed_expr.find(symbol_name, pos)) != std::string::npos) {
            parsed_expr.replace(pos, symbol_name.length(), symbol_rep);
            pos += symbol_rep.length();
        }
    }
    return SymEngine::parse(parsed_expr);
}

DenseMatrix compute_jacobian(const std::vector<ODE>& odes, const std::vector<RCP<const SymEngine::Symbol>>& state_vars) {
    size_t n = odes.size();
    DenseMatrix J(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            J.set(i,j,SymEngine::diff(odes[i].symbolic_expression.get_basic(), state_vars[j]));
        }
    }
    return J;
}

std::vector<RCP<const SymEngine::Symbol>> create_state_variable_symbols(const std::unordered_map<std::string, double>& state_var_values) {
    std::vector<RCP<const Symbol>> state_vars;
    for (const auto& var_name : state_var_values) {
        state_vars.push_back(SymEngine::rcp_static_cast<const Symbol>(symbol(var_name.first)));
    }
    return state_vars;
}

void identify_constants(ODE& ode, const std::unordered_set<std::string>& state_variables, std::unordered_set<std::string>& constants) {
    auto free_syms = SymEngine::free_symbols(*ode.symbolic_expression.get_basic());
    for (const auto& sym : free_syms) {
        std::string name = SymEngine::str(*sym);
        if (state_variables.find(name) == state_variables.end()) {
            constants.insert(name);
            ode.constants.insert(name);
        }
    }
}

std::vector<ODE> build_ode_system(const std::string& ode_str, std::unordered_set<std::string>& constants,
                                  symbol_table_t& symbol_table,
                                  std::unordered_map<std::string, double>& state_var_values,
                                  std::unordered_map<std::string, double>& constant_values) {
    std::vector<ODE> odes;
    std::unordered_set<std::string> state_variables;
    std::istringstream stream(ode_str);
    std::string line;
    std::unordered_map<std::string, RCP<const Basic>> symbols_map;
    SymEngine::Symbol t("t");
    // Extract state variables and parse ODEs
    while (std::getline(stream, line)) {
        std::istringstream line_stream(line);
        std::string name, eq, expr;
        line_stream >> name >> eq;
        std::getline(line_stream, expr);
        expr.erase(0, expr.find_first_not_of(" \t="));  // Remove leading spaces and '='

        // Replace all instances of '**' with '^' for compatibility with Exprtk
        replaceAll(expr, "**", "^");

        // Remove semicolons from the expression
        expr.erase(std::remove(expr.begin(), expr.end(), ';'), expr.end());

        std::string var_name = name.substr(1, name.find("_dt") - 1);  // Extract state variable name
        state_variables.insert(var_name);

        // Create symbolic variables
        symbols_map[var_name] = SymEngine::symbol(var_name);

        ODE ode = {name, expr, 0.0, {var_name}, {}};
        odes.push_back(ode);
    }

    // Identify constants and update state variables in each ODE
    for (auto& ode : odes) {
        std::unordered_set<std::string> vars;
        std::regex var_regex("[a-zA-Z_][a-zA-Z_0-9]*");
        auto words_begin = std::sregex_iterator(ode.expression_str.begin(), ode.expression_str.end(), var_regex);
        auto words_end = std::sregex_iterator();

        for (auto it = words_begin; it != words_end; ++it) {
            std::string var_name = it->str();
            if (state_variables.find(var_name) == state_variables.end() && var_name != "t") {  // Exclude state variables and known built-in variables
                vars.insert(var_name);
            }
        }

        //for (const auto& v : vars) {
        //    constants.insert(v);
        //    ode.constants.insert(v);
        //}

        // Update state variables for each ODE
        for (const auto& sv : state_variables) {
            if (ode.expression_str.find(sv) != std::string::npos) {
                ode.state_variables.push_back(sv);
            }
        }

        // Create symbolic representation of the expression
        ode.symbolic_expression = create_symbolic_expression(ode.expression_str, symbols_map);
        identify_constants(ode, state_variables, constants);
    }

    // Initialize ODE variables with initial values and add them to the symbol table
    for (const auto& var_name : state_variables) {
        state_var_values[var_name] = 0.0;  // Initial value can be set later
        symbol_table.add_variable(var_name, state_var_values[var_name]);
    }

    for (auto& ode : odes) {
        for (const auto& var_name : ode.state_variables) {
            ode.value = state_var_values[var_name];
        }
    }

    // Initialize constants with provided values and add them to the symbol table
    for (const auto& constant : constants) {
        auto it = constant_values.find(constant);
        if (it != constant_values.end()) {
            symbol_table.add_constant(constant, it->second);
        } else {
            symbol_table.add_constant(constant, 0.0);  // Default to 0.0 if not provided
        }
    }

    return odes;
}

std::function<void(const double, const Vector<double>&, Vector<double>&)> create_derivative_function(
        std::vector<ODE>& odes,
        symbol_table_t& symbol_table
) {
    std::vector<expression_t> expressions(odes.size());
    parser_t parser;

    for (size_t i = 0; i < odes.size(); ++i) {
        expressions[i].register_symbol_table(symbol_table);
        parser.compile(odes[i].expression_str, expressions[i]);
    }

    return [=](const double t, const Vector<double>& state, Vector<double>& derivatives) {
        // Update state variables in the symbol table
        for (size_t i = 0; i < state.size(); ++i) {
            symbol_table.get_variable(odes[i].state_variables[0])->ref() = state(i);
        }

        // Evaluate derivatives
        for (size_t i = 0; i < odes.size(); ++i) {
            derivatives(i) = expressions[i].value();
        }
    };
}

std::vector<std::string> extract_symbolic_expressions(const DenseMatrix& J) {
    std::vector<std::string> expressions;
    for (size_t i = 0; i < J.nrows(); ++i) {
        for (size_t j = 0; j < J.ncols(); ++j) {
            std::string expr_str = SymEngine::str(*J.get(i, j));
            replaceAll(expr_str, "**", "^");
            //std::cout << "Expression at (" << i << ", " << j << "): " << expr_str << std::endl;  // Debug output
            expressions.push_back(expr_str);
        }
    }
    return expressions;
}

std::function<void(const double, const Vector<double>&, Array<double>&)> create_jacobian_function(
        std::vector<ODE>& odes,
        const DenseMatrix& J,
        symbol_table_t &symbol_table) {

    std::vector<std::string> str_expressions = extract_symbolic_expressions(J);
    //std::cout << "Extracted expressions:" << std::endl;
    for (const auto& expr : str_expressions) {
        std::cout << expr << std::endl;
    }
    //std::cout << "Number of variables and constants: " << symbol_table.variable_count() << std::endl;
    //std::cout << "Get initial values:" << std::endl;
    for (const auto& var_name : odes[0].state_variables) {
        //std::cout << var_name << " = " << symbol_table.get_variable(var_name)->ref() << std::endl;
    }

    std::vector<expression_t> expressions(str_expressions.size());
    parser_t parser;

    for (size_t i = 0; i < str_expressions.size(); ++i) {
        expressions[i].register_symbol_table(symbol_table);
        if (!parser.compile(str_expressions[i], expressions[i])) {
            std::cerr << "Error: " << parser.error() << " in expression: " << str_expressions[i] << std::endl;
            //expressions[i].register_symbol_table(nullptr); // Ensure it is in a safe state
        } else {
            //std::cout << "Compiled expression: " << str_expressions[i] << std::endl;

            if (!is_valid(expressions[i])) {
                std::cerr << "Expression is not valid after compilation: " << str_expressions[i] << std::endl;
            }
        }
    }

    //std::cout << "Size of expressions vector: " << expressions.size() << std::endl;

    return [=](const double t, const Vector<double>& state, Array<double>& jacobian_values) {
        // Verify the state size matches expected size
        if (state.size() != odes.size()) {
            std::cerr << "State size does not match the number of ODEs" << std::endl;
            return;
        }

        // Ensure jacobian_values is correctly initialized
        //jacobian_values.resize(state.size(), std::vector<double>(state.size(), 0.0));
        jacobian_values.resize(state.size(), state.size());
        //std::cout << "jacobian_values initialized to size " << state.size() << "x" << state.size() << std::endl;

        //std::cout << "Constants:" << std::endl;
        for (const auto& constant : odes[0].constants) {
            auto var = symbol_table.get_variable(constant);
            if (var) {
                //std::cout << constant << " = " << var->ref() << std::endl;
            } else {
                std::cerr << "Constant " << constant << " not found in symbol table" << std::endl;
                return;
            }
        }

        for (size_t k = 0; k < state.size(); ++k) {
            auto var = symbol_table.get_variable(odes[k].state_variables[0]);
            symbol_table.get_variable(odes[k].state_variables[0])->ref() = state(k);
            if (var) {
                //var->ref() = state[k];
                //std::cout << "Updating variable " << odes[k].state_variables[0] << " with value " << state[k] << std::endl;
            } else {
                std::cerr << "Variable " << odes[k].state_variables[0] << " not found in symbol table" << std::endl;
                return;
            }

        }
        //std::cout << "Jacobian values:" << std::endl;
        for (size_t i = 0; i < state.size(); ++i) {
            for (size_t j = 0; j < state.size(); ++j) {
                // Update the symbol table with the current state for each evaluation
                /*
                for (size_t k = 0; k < state.size(); ++k) {
                    auto var = symbol_table.get_variable(odes[k].state_variables[0]);
                    symbol_table.get_variable(odes[k].state_variables[0])->ref() = state[k];
                    if (var) {
                        //var->ref() = state[k];
                        std::cout << "Updating variable " << odes[k].state_variables[0] << " with value " << state[k] << std::endl;
                    } else {
                        std::cerr << "Variable " << odes[k].state_variables[0] << " not found in symbol table" << std::endl;
                        return;
                    }
                }
                */

                size_t idx = i * state.size() + j;
                if (idx < expressions.size()) {
                    if (!is_valid(expressions[idx])) {
                        //std::cerr << "Expression at (" << i << ", " << j << ") is not valid." << std::endl;
                        continue;
                    }
                    try {
                        double value = expressions[idx].value();
                        //std::cout << "Evaluated expression at (" << i << ", " << j << "): " << value << std::endl;
                        jacobian_values(i, j) = value;
                    } catch (const std::exception& e) {
                        std::cerr << "Error evaluating expression at (" << i << ", " << j << "): " << e.what() << std::endl;
                    }
                } else {
                    std::cerr << "Index " << idx << " out of bounds for expressions." << std::endl;
                }
            }
        }

        // Output the final size of jacobian_values
        //std::cout << "Jacobian values size: " << jacobian_values.size() << std::endl;
    };
}

void getEigenValues(ODESystem& ode_system) {
    auto eigenvalues = SymEngine::eigen_values(ode_system.jacobian);
    std::cout << "Eigenvalues of the Jacobian matrix:" << std::endl;
    const SymEngine::FiniteSet& fs = static_cast<const SymEngine::FiniteSet&>(*eigenvalues);
    for (const auto& elem : fs.get_container()) {
        std::cout << *elem << std::endl;
    }
}

void printSparsity(const DenseMatrix& J) {
    std::cout << "Sparsity pattern of the Jacobian matrix:" << std::endl;
    for (size_t i = 0; i < J.nrows(); ++i) {
        for (size_t j = 0; j < J.ncols(); ++j) {
            if (J.get(i, j)->__eq__(*SymEngine::zero)) {
                std::cout << "0 ";
            } else {
                std::cout << "1 ";
            }
        }
        std::cout << std::endl;
    }
}

std::vector<std::pair<size_t, size_t>> getSparsityPattern(const DenseMatrix& J) {
    std::vector<std::pair<size_t, size_t>> sparsity;
    for (size_t i = 0; i < J.nrows(); ++i) {
        for (size_t j = 0; j < J.ncols(); ++j) {
            if (!J.get(i, j)->__eq__(*SymEngine::zero)) {
                sparsity.push_back(std::make_pair(i, j));
            }
        }
    }
    return sparsity;
}