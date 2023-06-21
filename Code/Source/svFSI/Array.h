
#include "Vector.h"
#include "utils.h"

#ifndef ARRAY_H 
#define ARRAY_H 

#include <algorithm>
#include <array>
#include <cstring>
#include <float.h>
#include <iostream>
#include <math.h>
#include <vector>

// If set then check Array indexes.
//
//#define Array_check_enabled

//-------
// Array 
//-------
// The Array template class implements a simple interface to 2D arrays.
//
template<typename T>
class Array 
{
  public:
    // Some variables used to monitor Array object allocation.
    static int id;
    static int num_allocated;
    static int active;
    static double memory_in_use;
    static double memory_returned;
    static void memory(const std::string& prefix="");
    static void stats(const std::string& prefix="");
    static bool write_enabled;

    Array() 
    {
      id += 1;
      nrows_ = 0;
      ncols_ = 0;
      size_ = 0;
      data_ = nullptr;
      num_allocated += 1;
      active += 1;
    };

    Array(const int num_rows, const int num_cols)
    {
      // [NOTE] This is tfu but need to mimic fortran.
      nrows_ = num_rows;
      ncols_ = num_cols;

      if ((num_rows <= 0) || (num_cols <= 0)) { 
        return;
      }

      allocate(num_rows, num_cols);
      num_allocated += 1;
      active += 1;
    }

    Array(const int num_rows, const int num_cols, T* data)
    {
      data_reference_ = true;
      nrows_ = num_rows;
      ncols_ = num_cols;
      size_ = nrows_ * ncols_;
      data_ = data;
    }

    Array(std::initializer_list<std::initializer_list<T>> lst)
    {
      int num_rows = 0;
      int num_cols = 0;

      for (auto &row : lst) {
        num_rows += 1;
        if (num_cols == 0) {
          num_cols = row.size();
        } else if (row.size() != num_cols) {
          throw std::runtime_error("[Array initializer_list ctor] Column data is of unequal size.");
        }
      }

      allocate(num_rows, num_cols);

      int row = 0;

      for (auto &row_data : lst) {
        auto it = row_data.begin();
        for (int col = 0; col < row_data.size(); col++, it++) {
          data_[row + col*nrows_] = *it;
        }
        row += 1;
      }

      num_allocated += 1;
      active += 1;
    }

    //-------------
    // Array copy
    //-------------
    //
    Array(const Array &rhs)
    {
      if ((rhs.nrows_ <= 0) || (rhs.ncols_ <= 0)) { 
        return;
      }

      allocate(rhs.nrows_, rhs.ncols_);

      memcpy(data_, rhs.data_, size_*sizeof(T));
      num_allocated += 1;
      active += 1;
    }

    // Array assignment 
    //
    Array& operator=(const Array& rhs)
    {
      if (this == &rhs) {
        return *this;
      }

      if ((rhs.nrows_ <= 0) || (rhs.ncols_ <= 0)) { 
        return *this;
      }

      if ((nrows_ != rhs.nrows_) || (ncols_ != rhs.ncols_)) {     
        clear();
        allocate(rhs.nrows_, rhs.ncols_);
      } 

      memcpy(data_, rhs.data_, sizeof(T) * size_);
      return *this;
    }

    Array& operator=(const double value)
    {
      for (int i = 0; i < size_; i++) {
        data_[i] = value;
      }
      return *this;
    }

    ~Array() 
    {
      if (data_ != nullptr) {
        #if Array_gather_stats
        memory_in_use -= sizeof(T)*size_;
        memory_returned += sizeof(T)*size_;
        active -= 1;
        #endif
        if (!data_reference_) {
          delete [] data_;
        }
        data_ = nullptr;
        size_ = 0;
        id -= 1;
        nrows_ = 0;
        ncols_ = 0;
       }
    };

    friend std::ostream& operator << (std::ostream& out, const Array<T>& lhs)
    {
      for (int i = 0; i < lhs.size(); i++) {
        out << lhs.data_[i];
        if (i != lhs.size()-1) {
          out << ", ";
         }
       }
       return out;
    }

    //-------
    // clear 
    //-------
    // Free the array data. 
    //
    // This is to replicate the Fortran DEALLOCATE().
    //
    void clear()
    {
      if (data_ != nullptr) {
        if (data_reference_) {
          throw std::runtime_error("[Array] Can't clear an Array with reference data.");
        }
        delete [] data_;
        #if Array_gather_stats
        memory_in_use -= sizeof(T) * size_;;
        memory_returned += sizeof(T) * size_;;
        #endif
      }
      nrows_ = 0;
      ncols_ = 0;
      size_ = 0;
      data_ = nullptr;
    }

    //----------
    // ncols
    //----------
    int ncols() const
    {
      return ncols_;
    }

    //----------
    // nrows
    //----------
    int nrows() const
    {
      return nrows_;
    }

    //-------
    // Print
    //-------
    // Print by rows and columns.
    //
    void print(const std::string& label)
    {
      printf("%s %d x %d \n", label.c_str(), nrows_, ncols_);
      for (int i = 0; i < nrows_; i++) {
        printf("[ ");
        for (int j = 0; j < ncols_; j++) {
          printf("%.16e", value(i,j));
          if (j == ncols_-1) {
            printf(" ]");
          } else {
            printf(", ");
          }
        }
        printf("\n");
      }
    }

    //--------
    // resize
    //--------
    // Resize the array.
    //
    void resize(const int num_rows, const int num_cols) 
    {
      // [NOTE] This is tfu but need to mimic fortran.
      nrows_ = num_rows;
      ncols_ = num_cols;

      if ((num_rows <= 0) || (num_cols <= 0)) { 
        return;
      }

      if (data_ != nullptr) {
        if (data_reference_) {
          throw std::runtime_error("[Array] Can't resize an Array with reference data.");
        }
        delete [] data_;
        data_ = nullptr;
        size_ = 0;
        nrows_ = 0;
        ncols_ = 0;
        #if Array_gather_stats
        memory_in_use -= sizeof(T) * size_;;
        memory_returned += sizeof(T) * size_;;
        #endif
      }

      allocate(num_rows, num_cols);
    }

    //------
    // size
    //------
    int size() const
    {
      return size_;
    }

    int msize() const
    {
      return size_ * sizeof(T);
    }

    // Test if an array has rows or columns set
    // but no data.
    //
    // [NOTE] This is tfu but need to mimic fortran.
    //
    bool allocated()
    {
      if ((nrows_ > 0) || (ncols_ > 0)) { 
        return true;
      }

      return false;
    }

    //------
    // read
    //------
    void read(const std::string& file_name) 
    {
      auto fp = fopen(file_name.c_str(), "rb");
      int size;
      fread(&size, sizeof(int), 1, fp);
      fread(data_, sizeof(double), size, fp);
      fclose(fp);
    }

    //-------
    // write
    //-------
    //
    void write(const std::string& label, bool memory=true, T offset={}) const
    {
      if (!write_enabled) {
        return;
      }

      auto file_prefix = build_file_prefix(label);
      auto file_name = file_prefix + "_cm.bin";

      // Write binary file.
      //
      auto fp = fopen(file_name.c_str(), "wb");
      fwrite(&size_, sizeof(int), 1, fp);
      fwrite(data_, sizeof(T), size_, fp);
      fclose(fp);
    }

    //--------
    // append
    //--------
    // Append data to a file.
    //
    void append(const std::string& label, bool memory=true, T offset={}) const
    {
      if (!write_enabled) {
        return;
      }

      auto file_prefix = build_file_prefix(label);
      auto file_name = file_prefix + "_cm.bin";

      // Append binary file.
      //
      auto fp = fopen(file_name.c_str(), "ab");
      fwrite(data_, sizeof(T), size_, fp);
      fclose(fp);
    }

    const T& operator()(const int i) const
    {
      if (i >= size_) {
        throw std::runtime_error("[Array(i)] Index " + std::to_string(i) + " is out of bounds.");
      }
      return data_[i];
    }

    T& operator()(const int i)
    {
      if (i >= size_) {
        throw std::runtime_error("[Array(i)] Index " + std::to_string(i) + " is out of bounds.");
      }
      return data_[i];
    }

    //-----------------
    // (i,j) operators
    //-----------------
    //
    const T& operator()(const int row, const int col) const
    {
      #ifdef Array_check_enabled
      check_index(row, col);
      #endif
      return data_[row + col*nrows_];
    }

    T& operator()(const int row, const int col)
    {
      #ifdef Array_check_enabled
      check_index(row, col);
      #endif
      return data_[row + col*nrows_];
    }

    //------
    // col 
    //------
    // Get a column from the array as a Vector.
    //
    Vector<T> col(const int col, const std::array<int,2>& range={-1,-1}) const
    {
      #ifdef Array_check_enabled
      check_index(0, col);
      #endif
      int start_row, end_row;

      if (range[0] != -1) {
        start_row = range[0];
        end_row = range[1];
        #ifdef Array_check_enabled
        check_index(start_row, col);
        check_index(end_row, col);
        #endif
      } else {
        start_row = 0;
        end_row = nrows_;
      }

      Vector<T> values(nrows_);
      for (int row = start_row; row  < end_row; row++) {
        values[row] = value(row,col);
      }
      return values;
    }

    //------
    // rcol
    //------
    Vector<T> rcol(const int col) const
    {
      Vector<T> vector_col(nrows_, &data_[col*nrows_]);
      return vector_col;
    }

    //----------
    // col_data
    //----------
    // Return a pointer to the internal data for the given column.
    //
    T* col_data(const int col) { 
      return &data_[col*nrows_];
    }

    //------
    // cols
    //------
    // Get one or more columns from the array as Vectors.
    //
    Vector<Vector<T>> cols(const Vector<int>& columns) const
    {
      Vector<Vector<T>> values(columns.size());

      for (int i = 0; i < columns.size(); i++) {
        int j = columns[i];
        #ifdef Array_check_enabled
        check_index(0, j);
        #endif
        values[i] = col(j); 
      }

      return values;
    }

    //-----
    // row 
    //-----
    // Get a row from the array as a Vector.
    //
    Vector<T> row(const int row, const std::array<int,2>& range={-1,-1})
    {
      #ifdef Array_check_enabled
      check_index(row, 0);
      #endif

      Vector<T> values(ncols_);

      for (int col = 0; col < ncols_; col++) {
        values[col] = value(row,col);
      }

      return values;
    }

    //------
    // rows
    //------
    // Get row values from the array using an index of columns.
    //
    Vector<T> rows(const int row, const Vector<int>& indexes) const
    {
      #ifdef Array_check_enabled
      check_index(row, 0);
      #endif

      Vector<T> values(indexes.size());

      for (int i = 0; i < indexes.size(); i++) {
        int col = indexes[i];
        #ifdef Array_check_enabled
        check_index(row, col);
        #endif
        values[i] = value(row,col);
      }

      return values;
    }

    //------
    // rows
    //------
    // Get a set of rows using a start and end row index.
    //
    Array<T> rows(const int start_row, const int end_row) const
    {
      #ifdef Array_check_enabled
      check_index(start_row, 0);
      check_index(end_row, 0);
      #endif
      int num_rows = end_row - start_row + 1;
      Array<T> values(num_rows, ncols_ );

      for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < ncols_; col++) {
          values(row,col) = data_[(row+start_row) + col*nrows_];
        }
      }
      return values;
    }

    //------
    // rows
    //------
    // Get a range of rows for the given column.
    //
    Vector<T> rows(const int start_row, const int end_row, const int col) const
    {
      #ifdef Array_check_enabled
      check_index(start_row, 0);
      check_index(end_row, 0);
      #endif
      int num_rows = end_row - start_row + 1;
      Vector<T> values(num_rows);

      for (int row = 0; row < num_rows; row++) {
        values(row) = data_[(row+start_row) + col*nrows_];
      }
      return values;
    }

    //-------
    // values
    //--------
    // Return a vector of values for a range of rows and columns.
    //
    std::vector<T> values(const std::array<int,2>& rows, const std::array<int,2>& cols, const int stride=1) const
    {
      std::vector<T> values;

      for (int i = rows[0]; i <= rows[1]; i += stride) {
        for (int j = cols[0]; j <= cols[1]; j += stride) {
          values.push_back(value(i,j));
        }
      }
    
      return values;
    }

    //---------
    // set_col
    //---------
    // Set column values from a vector of values.
    //
    void set_col(const int col, const Vector<T>& values, const std::array<int,2>& range={-1,-1}) 
    {
      #ifdef Array_check_enabled
      check_index(0, col);
      #endif
      int start_row, end_row;

      if (range[0] != -1) {
        start_row = range[0];
        end_row = range[1];
        #ifdef Array_check_enabled
        check_index(start_row, col);
        check_index(end_row, col);
        #endif
      } else {
        start_row = 0;
        end_row = nrows_;
      }

      for (int row = 0; row < values.size(); row++) {
        data_[row + col*nrows_] = values[row];
      }
    }

    //---------
    // set_row
    //---------
    // Set row values from a vector of values.
    //
    void set_row(const int row, const Vector<T>& values) const
    {
      #ifdef Array_check_enabled
      check_index(row, 0);
      #endif
      for (int col = 0; col < values.size(); col++) {
        data_[row + col*nrows_] = values[col];
      }
    }

    //---------
    // set_row
    //---------
    // Set row values from an initializer list {}..
    //
    void set_row(const int row, std::initializer_list<T> row_data) const
    {
      auto it = row_data.begin();
      for (int col = 0; col < ncols_; col++, it++) {
          data_[row + col*nrows_] = *it;
      }
    }

    //---------
    // set_row
    //---------
    // Set row values from a scalar. 
    //
    void set_row(const int row, const T value) const
    {
      #ifdef Array_check_enabled
      check_index(row, 0);
      #endif
      for (int col = 0; col < ncols_; col++) {
        data_[row + col*nrows_] = value;
      }
    }

    //---------
    // set_row
    //---------
    // Set row values from a vector of values at the given column offset.
    //
    void set_row(const int row, int offset, const Vector<T>& values) const
    {
      #ifdef Array_check_enabled
      check_index(row, 0);
      check_index(row, offset+values.size()-1);
      #endif
      for (int col = 0; col < values.size(); col++) {
        int ocol = col + offset;
        #ifdef Array_check_enabled
        check_index(row, ocol);
        #endif
        data_[row + ocol*nrows_] = values[col];
      }
    }

    //----------
    // set_rows
    //----------
    // Set a set of rows using a start and end row index.
    //
    void set_rows(const int start_row, const int end_row, const Array<T>& values) const
    { 
      #ifdef Array_check_enabled
      check_index(start_row, 0);
      check_index(end_row, 0); 
      #endif
      int num_rows = end_row - start_row + 1;

      if (ncols_ != values.ncols_) {
        throw std::runtime_error("[Array set_rows] Arrays have different column sizes. ");
      }
      
      for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < ncols_; col++) {
          data_[(row+start_row) + col*nrows_] = values(row,col);;
        }
      }
    }

    ///////////////////////////////////////////////////
    //  M a t h e m a t i c a l   o p e r a t o r s  //
    ///////////////////////////////////////////////////

    //-----
    // max
    //-----
    T max() const
    {
      T max_v = data_[0];
      for (int i = 1; i < size_; i++) {
        if (data_[i] > max_v) {
          max_v = data_[i];
        }
      }
      return max_v;
    }

    //-----
    // min
    //-----
    T min() const
    {
      T min_v = data_[0];
      for (int i = 1; i < size_; i++) {
        if (data_[i] < min_v) {
          min_v = data_[i];
        }
      }
      return min_v;
    }

    // Add arrays.
    //
    Array<T> operator+(const Array<T>& array) const
    {
      if ((nrows_ != array.nrows_) || (ncols_ != array.ncols_)) {
        throw std::runtime_error("[Array addition] Arrays have diffent shapes. ");
      }
      Array<T> result(nrows_, ncols_);
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          result(i, j) = data_[i + j*nrows_] + array(i,j);
        }
      }
      return result;
    }

    // Subtract arrays.
    //
    Array<T> operator-(const Array<T>& array) const
    {
      if ((nrows_ != array.nrows_) || (ncols_ != array.ncols_)) {
        throw std::runtime_error("[Array subtraction]  Arrays have diffent shapes. ");
      }
      Array<T> result(nrows_, ncols_);
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          result(i, j) = data_[i + j*nrows_] - array(i,j);
        }
      }
      return result;
    }

    //-----------
    // C = A * B
    //-----------
    // Multiply two arrays component-wise to reproduce Fortran.
    //
    Array<T> operator*(const Array<T>& array) const
    {
      if ((nrows_ != array.nrows_) || (ncols_ != array.ncols_)) {
        throw std::runtime_error("[Array multiply] Arrays have diffent shapes. ");
      }

      Array<T> result(nrows_, ncols_);

      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          result(i, j) = (*this)(i,j) * array(i,j);
        }
      }

      return result;
    }

    //-----------
    // C = A / B
    //-----------
    // Divide one array by another component-wise.
    //
    Array<T> operator / (const Array<T>& array) const
    {
      Array<T> result(nrows_, array.ncols_);
      if ((nrows_ != array.nrows_) || (ncols_ != array.ncols_)) {
        throw std::runtime_error("[Array divide] Arrays number of columns or number of rows are not equal.");
      }

      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          result(i, j) = (*this)(i,j) / array(i,j);
        }
      }
      return result;
    }

    // Compound add assignment. 
    //
    Array<T> operator+=(const Array<T>& array) const
    {
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          data_[i + j*nrows_] += array(i,j);
        }
      }
      return *this;
    }

    // Compound subtract assignment. 
    //
    Array<T> operator-=(const Array<T>& array) const
    {
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          data_[i + j*nrows_] -= array(i,j);
        }
      }
      return *this;
    }

    // Compound multiply assignment. 
    //
    Array<T> operator*=(const Array<T>& array) const
    {
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          data_[i + j*nrows_] *= array(i,j);
        }
      }
      return *this;
    }

    // Multiply by a scalar.
    //
    Array<T> operator*(const T value) const
    {
      Array<T> result(nrows_, ncols_);
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          result(i,j) = value * data_[i + j*nrows_];
        }
      }
      return result;
    }

    friend const Array<T> operator*(const T value, const Array& rhs)
    {
      Array<T> result(rhs.nrows_, rhs.ncols_);
      for (int j = 0; j < rhs.ncols_; j++) {
        for (int i = 0; i < rhs.nrows_; i++) {
          result(i,j) = value * rhs.data_[i + j*rhs.nrows_];
        }
      }
      return result;
    }

    // A / s
    //
    // Divide by a scalar.
    //
    Array<T> operator / (const T value) const
    {
      if (value == 0.0) {
        throw std::runtime_error(+"Array Divide by zero.");
      }
      Array<T> result(nrows_, ncols_);
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          result(i,j) = data_[i + j*nrows_] / value;
        }
      }
      return result;
    }

    //-------
    // s / A 
    //-------
    // Divide scalar by array.
    //
    friend const Array<T> operator / (const T value, const Array& rhs)
    {
      Array<T> result(rhs.nrows_, rhs.ncols_);
      for (int j = 0; j < rhs.ncols_; j++) {
        for (int i = 0; i < rhs.nrows_; i++) {
          result(i,j) = value / rhs.data_[i + j*rhs.nrows_];
        }
      }
      return result;
    }

    //-------
    // A - s
    //-------
    // Subtract a scalar.
    //
    Array<T> operator-(const T value) const
    { 
      Array<T> result(nrows_, ncols_); 
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) { 
          result(i,j) = data_[i + j*nrows_] - value;
        }
      }
      return result;
    }

    // Compound add assignment. 
    //
    Array<T> operator+=(const T value) const
    {
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          data_[i + j*nrows_] += value;
        }
      }
      return *this;
    }

    //--------
    // negate 
    //--------
    //
    Array<T> operator-() const 
    {
      Array<T> result(nrows_, ncols_);
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          result(i,j) = -(data_[i + j*nrows_]);
        }
      }
      return result;
    }

    // Compound subtract assignment. 
    //
    Array<T> operator-=(const T value) const
    {
      for (int j = 0; j < ncols_; j++) {
        for (int i = 0; i < nrows_; i++) {
          data_[i + j*nrows_] -= value;
        }
      }
      return *this;
    }

    //-------
    // s - A
    //-------
    // Subtract a scalar.
    //
    friend const Array<T> operator-(const T value, const Array& rhs)
    { 
      Array<T> result(rhs.nrows_, rhs.ncols_);
      for (int j = 0; j < rhs.ncols_; j++) {
        for (int i = 0; i < rhs.nrows_; i++) { 
          result(i,j) = value - rhs.data_[i + j*rhs.nrows_];
        }
      }
      return result;
    }

    //------
    // abs
    //------
    //
    friend Array<T> abs(const Array& rhs)
    {
      Array<T> result(rhs.nrows_, rhs.ncols_);
      for (int j = 0; j < rhs.ncols_; j++) {
        for (int i = 0; i < rhs.nrows_; i++) {
          result(i,j) = fabs(rhs.data_[i + j*rhs.nrows_]);
        }
      }
      return result;
    }

    //-----
    // max
    //-----
    friend T max(const Array& arg) 
    {
      T max_v = arg.data_[0];
      for (int i = 1; i < arg.size_; i++) {
        if (arg.data_[i] > max_v) {
          max_v = arg.data_[i];
        }
      }
      return max_v;
    }

    //------
    // sqrt
    //------
    //
    friend Array<T> sqrt(const Array& arg)
    {
      Array<T> result(arg.nrows_, arg.ncols_);
      for (int j = 0; j < arg.ncols_; j++) {
        for (int i = 0; i < arg.nrows_; i++) {
          result(i,j) = sqrt(arg.data_[i + j*arg.nrows_]);
        }
      }
      return result;
    }

    //---------
    // sum_row
    //---------
    // Compute the sum of a row of valuesr
    //
    T sum_row(const int row) const
    {
      #ifdef Array_check_enabled
      check_index(row, 0);
      #endif
      T sum {};
      for (int col = 0; col < ncols_; col++) {
        sum += data_[row + col*nrows_];
      }
      return sum;
    }

    T* data() const { 
      return data_;
    }

  private:

    //----------
    // allocate
    //----------
    // Allocate memory for array data.
    //
    void allocate(const int num_rows, const int num_cols)
    {
      data_reference_ = false;
      nrows_ = num_rows;
      ncols_ = num_cols;
      size_ = nrows_ * ncols_;
      #if Array_gather_stats
      memory_in_use += sizeof(T)*size_;
      #endif

      if (data_ != nullptr) { 
        //throw std::runtime_error(+"[Array] Allocating when data is not nullptr.");
      }

      if (size_ != 0) {
        data_ = new T [size_];
        memset(data_, 0, sizeof(T)*size_);
      }
    }

    //-------
    // value
    //-------
    // Get a value from data_[].
    //
    inline T value(const int row, const int col) const
    {
      return data_[row + col*nrows_];
    }

    //-------------
    // check_index
    //-------------
    //
    void check_index(const int row, const int col) const
    { 
      if (data_ == nullptr) { 
        //throw std::runtime_error(+"Accessing null data in Array.");
        return;
      }
      if ((row < 0) || (row >= nrows_) || (col < 0) || (col >= ncols_)) {
        auto nr_str = std::to_string(nrows_);
        auto nc_str = std::to_string(ncols_);
        auto dims = nr_str + " x " + nc_str;
        auto index_str = " " + std::to_string(row) + "," + std::to_string(col) + " ";
        throw std::runtime_error(+"Index (row,col)=" + index_str + " is out of bounds for " + 
            dims + " array.");
      }
    }

    //-------------
    // Member data
    //-------------
    //
    int nrows_ = 0;
    int ncols_ = 0;
    int size_ = 0;
    bool data_reference_ = false;
    T *data_ = nullptr;
};

#endif

