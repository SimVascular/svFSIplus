import assert from 'node:assert/strict';
import test from 'node:test';
import { evaluate, GitHub, reconcile } from './review-policy.mjs';

const author = { id: 1, login: 'author', type: 'User' };
const specialist = { id: 2, login: 'specialist', type: 'User' };
const maintainer = { id: 3, login: 'maintainer', type: 'User' };
const outsider = { id: 4, login: 'outsider', type: 'User' };
const label = { id: 9482575472, name: 'OOP Refactor' };
const team = new Set([author.id, specialist.id]);
const writers = new Set([author.id, specialist.id, maintainer.id]);

function review(id, user, state = 'APPROVED') {
  return { id, user, state, submitted_at: new Date(Date.UTC(2026, 9, 12, 10, id)).toISOString() };
}

function inputs(overrides = {}) {
  return {
    pr: { user: author, labels: [] },
    reviews: [],
    timeline: [],
    policyChanged: false,
    ...overrides,
  };
}

function oop(overrides = {}) {
  return inputs({
    pr: { user: author, labels: [label] },
    timeline: [{ id: 10, event: 'labeled', label, actor: author }],
    ...overrides,
  });
}

test('ordinary PRs need two distinct reviewers with write permission', () => {
  assert.equal(evaluate(inputs({ reviews: [review(1, specialist)] }), team, writers).passed, false);
  assert.equal(evaluate(inputs({ reviews: [review(1, specialist), review(2, maintainer)] }), team, writers).passed, true);
  assert.equal(evaluate(inputs({ reviews: [review(1, specialist), review(2, specialist)] }), team, writers).passed, false);
  assert.equal(evaluate(inputs({ reviews: [review(1, specialist), review(2, outsider)] }), team, writers).passed, false);
});

test('an authorized OOP label permits one independent team approval', () => {
  const result = evaluate(oop({ reviews: [review(1, specialist)] }), team, writers);
  assert.equal(result.passed, true);
  assert.equal(result.required, 1);
  assert.equal(result.removeLabel, false);
  assert.equal(evaluate(oop({ reviews: [review(1, author)] }), team, writers).passed, false);
  assert.equal(evaluate(oop({ reviews: [review(1, maintainer)] }), team, writers).passed, false);
  assert.equal(evaluate(oop({ reviews: [review(1, specialist)] }), team, new Set()).passed, false);
});

test('the labeler can approve another author’s PR', () => {
  const result = evaluate(oop({
    pr: { user: maintainer, labels: [label] },
    reviews: [review(1, author)],
  }), team, writers);
  assert.equal(result.passed, true);
});

test('unauthorized or untraceable labels are removed without granting an exception', () => {
  for (const actor of [outsider, null]) {
    const result = evaluate(oop({
      timeline: [{ id: 10, event: 'labeled', label, actor }],
      reviews: [review(1, specialist)],
    }), team, writers);
    assert.equal(result.removeLabel, true);
    assert.equal(result.required, 2);
    assert.equal(result.passed, false);
  }
  assert.equal(evaluate(oop({ timeline: [] }), team, writers).removeLabel, true);
});

test('reapplying a label cannot reuse an earlier authorized actor', () => {
  const result = evaluate(oop({
    timeline: [
      { id: 10, event: 'labeled', label, actor: author },
      { id: 11, event: 'unlabeled', label, actor: outsider },
      { id: 12, event: 'labeled', label, actor: outsider },
    ],
    reviews: [review(1, specialist)],
  }), team, writers);
  assert.equal(result.removeLabel, true);
  assert.equal(result.passed, false);
});

test('label removal and current team membership restore the ordinary policy', () => {
  assert.equal(evaluate(inputs({ reviews: [review(1, specialist)] }), team, writers).required, 2);
  assert.equal(evaluate(oop({ reviews: [review(1, specialist)] }), new Set([specialist.id]), writers).required, 2);
  assert.equal(evaluate(oop({ reviews: [review(1, specialist)] }), new Set([author.id]), writers).passed, false);
});

test('a recreated label with the same name cannot acquire the exception', () => {
  const result = evaluate(oop({ pr: { user: author, labels: [{ ...label, id: 99 }] } }), team, writers);
  assert.equal(result.required, 2);
  assert.equal(result.removeLabel, true);
});

test('comments preserve approvals while dismissals and change requests replace them', () => {
  for (const [state, passed] of [['COMMENTED', true], ['PENDING', true], ['DISMISSED', false], ['CHANGES_REQUESTED', false]]) {
    const result = evaluate(oop({ reviews: [review(1, specialist), review(2, specialist, state)] }), team, writers);
    assert.equal(result.passed, passed, state);
  }
  assert.equal(evaluate(oop({ reviews: [review(1, specialist, 'DISMISSED'), review(2, specialist)] }), team, writers).passed, true);
});

test('an eligible change request blocks approval even from a different reviewer', () => {
  const reviews = [review(1, specialist), review(2, maintainer, 'CHANGES_REQUESTED')];
  assert.equal(evaluate(oop({ reviews }), team, writers).passed, false);
  assert.equal(evaluate(oop({ reviews }), team, new Set([specialist.id])).passed, true);
});

test('review submission order takes precedence over draft creation order', () => {
  const reviews = [review(20, specialist), {
    ...review(10, specialist, 'CHANGES_REQUESTED'), submitted_at: '2026-10-12T11:00:00Z',
  }];
  assert.equal(evaluate(oop({ reviews }), team, writers).passed, false);
});

test('bots and deleted reviewers never supply an approval', () => {
  for (const user of [null, { ...specialist, type: 'Bot' }]) {
    assert.equal(evaluate(oop({ reviews: [review(1, user)] }), team, writers).passed, false);
  }
});

test('changes to policy enforcement require two approvals even with the label', () => {
  const result = evaluate(oop({ policyChanged: true, reviews: [review(1, specialist)] }), team, writers);
  assert.equal(result.required, 2);
  assert.equal(result.passed, false);
  assert.equal(evaluate(oop({ policyChanged: true, reviews: [review(1, specialist), review(2, maintainer)] }), team, writers).passed, true);
});

const repoPath = '/repos/SimVascular/svMultiPhysics';
const teamPath = '/orgs/SimVascular/teams/refactoring-team';
const head = 'a'.repeat(40);

function pull(number = 645, overrides = {}) {
  return {
    number, state: 'open', user: author, labels: [], draft: false,
    head: { sha: head, repo: { full_name: 'contributor/svMultiPhysics' } },
    base: { ref: 'main', sha: 'b'.repeat(40) }, updated_at: '2026-10-12T10:00:00Z',
    ...overrides,
  };
}

// Model the HTTP boundary; policy evaluation and API pagination run unchanged.
function server(prs = [pull()]) {
  const state = {
    prs: structuredClone(prs), members: [author, specialist], teamId: 20019155,
    reviews: new Map(prs.map(pr => [pr.number, [review(1, specialist), review(2, maintainer)]])),
    timelines: new Map(prs.map(pr => [pr.number, oop().timeline])),
    files: new Map(prs.map(pr => [pr.number, []])), checks: [], writes: [],
    permissions: new Map([[author.id, 'write'], [specialist.id, 'write'], [maintainer.id, 'admin'], [outsider.id, 'read']]),
    before: () => {},
  };
  const fetcher = async (input, options) => {
    const url = new URL(input);
    const path = url.pathname;
    const method = options.method;
    const body = options.body ? JSON.parse(options.body) : undefined;
    if (method !== 'GET') state.writes.push({ path, method, body });
    const override = state.before({ path, method, url, body });
    if (override) return override;
    const page = array => array.slice((Number(url.searchParams.get('page') ?? 1) - 1) * 100,
      Number(url.searchParams.get('page') ?? 1) * 100);
    let data;
    if (path === `${repoPath}/pulls`) {
      data = page(state.prs.filter(pr => pr.state === 'open' && pr.base.ref === 'main'));
    } else if (path === teamPath) {
      data = { id: state.teamId, slug: 'refactoring-team' };
    } else if (path === `${teamPath}/members`) {
      data = page(state.members);
    } else if (path.includes('/collaborators/')) {
      const login = path.split('/').at(-2);
      const user = [author, specialist, maintainer, outsider].find(item => item.login === login);
      data = { permission: state.permissions.get(user.id), user };
    } else if (path.includes('/commits/')) {
      const sha = path.split('/').at(-2);
      data = { check_runs: page(state.checks.filter(check => check.head_sha === sha)) };
    } else if (path === `${repoPath}/check-runs` && method === 'POST') {
      data = { ...body, id: 1000 + state.checks.length, app: { id: 5288861 } };
      state.checks.push(data);
    } else if (path.startsWith(`${repoPath}/check-runs/`) && method === 'PATCH') {
      data = state.checks.find(check => check.id === Number(path.split('/').at(-1)));
      // Exercise the policy without relying on any mutation of completed checks.
      if (data.status === 'completed') return Response.json({}, { status: 422 });
      Object.assign(data, body);
      if (body.status === 'in_progress') delete data.conclusion;
    } else {
      const match = path.match(/\/(pulls|issues)\/(\d+)(.*)$/);
      assert.ok(match, `Unhandled request: ${method} ${path}`);
      const number = Number(match[2]);
      const pr = state.prs.find(item => item.number === number);
      switch (match[3]) {
        case '': data = pr; break;
        case '/reviews': data = page(state.reviews.get(number)); break;
        case '/timeline': data = page(state.timelines.get(number)); break;
        case '/files': data = page(state.files.get(number)); break;
        case '/labels/OOP%20Refactor':
          assert.equal(method, 'DELETE');
          pr.labels = pr.labels.filter(item => item.name !== 'OOP Refactor');
          data = pr.labels;
          break;
        default: assert.fail(`Unhandled request: ${method} ${path}`);
      }
    }
    return Response.json(data);
  };
  return { state, api: new GitHub('test-token', fetcher) };
}

test('dry-run evaluates fork PRs without creating checks or removing labels', async () => {
  const { api, state } = server([pull(645, { labels: [label] })]);
  state.timelines.set(645, [{ id: 10, event: 'labeled', label, actor: outsider }]);
  const results = await reconcile(api);
  assert.equal(results[0].removeLabel, true);
  assert.equal(results[0].required, 2);
  assert.deepEqual(state.writes, []);
  assert.deepEqual(state.prs[0].labels, [label]);
});

test('the app publishes the check on the fork head rather than main', async () => {
  const { api, state } = server([pull(645, { labels: [label] })]);
  state.reviews.set(645, [review(1, specialist)]);
  await reconcile(api, { apply: true });
  assert.equal(state.checks.length, 1);
  assert.equal(state.checks[0].head_sha, head);
  assert.equal(state.checks[0].name, 'review-policy');
  assert.equal(state.checks[0].conclusion, 'success');
});

test('shared commit checks pass only when every associated PR satisfies its policy', async () => {
  const { api, state } = server([pull(645, { labels: [label] }), pull(646)]);
  state.reviews.set(645, [review(1, specialist)]);
  state.reviews.set(646, [review(1, specialist)]);
  await reconcile(api, { apply: true });
  assert.equal(state.checks.length, 1);
  assert.equal(state.checks[0].conclusion, 'failure');
  assert.match(state.checks[0].output.summary, /645/);
  assert.match(state.checks[0].output.summary, /646/);
});

test('apply removes unauthorized labels and evaluates the ordinary policy', async () => {
  const { api, state } = server([pull(645, { labels: [label] })]);
  state.timelines.set(645, [{ id: 10, event: 'labeled', label, actor: outsider }]);
  state.reviews.set(645, [review(1, specialist)]);
  await reconcile(api, { apply: true });
  assert.deepEqual(state.prs[0].labels, []);
  assert.equal(state.checks[0].conclusion, 'failure');
});

test('policy changes cannot grant themselves a one-review exception', async () => {
  const { api, state } = server([pull(645, { labels: [label] })]);
  state.reviews.set(645, [review(1, specialist)]);
  state.files.set(645, [{ filename: '.github/workflows/review-policy-events.yml', status: 'modified' }]);
  await reconcile(api, { apply: true });
  assert.equal(state.checks[0].conclusion, 'failure');
});

test('any workflow change requires two approvals because workflows can request the environment', async () => {
  for (const filename of ['.github/workflows/new.yml', '.github/workflows/tests.yml']) {
    const { api, state } = server([pull(645, { labels: [label] })]);
    state.reviews.set(645, [review(1, specialist)]);
    state.files.set(645, [{ filename, status: 'modified' }]);
    await reconcile(api, { apply: true });
    assert.equal(state.checks[0].conclusion, 'failure', filename);
  }
});

test('API errors replace previously successful checks with failure', async () => {
  const { api, state } = server();
  state.checks.push({ id: 42, head_sha: head, name: 'review-policy', app: { id: 5288861 }, status: 'completed', conclusion: 'success' });
  state.before = ({ path }) => path === `${teamPath}/members` ? Response.json({}, { status: 403 }) : undefined;
  await assert.rejects(reconcile(api, { apply: true }), /403/);
  assert.equal(state.checks[0].conclusion, 'success');
  assert.equal(state.checks.at(-1).conclusion, 'failure');
  assert.equal(state.checks.length, 2);
  assert.equal(state.writes.some(write => write.path.endsWith('/42')), false);
});

test('unchanged completed results need no writes on scheduled refreshes', async () => {
  for (const reviews of [[], [review(1, specialist), review(2, maintainer)]]) {
    const { api, state } = server();
    state.reviews.set(645, reviews);
    await reconcile(api, { apply: true, detailsUrl: 'https://example.org/runs/1' });
    const original = structuredClone(state.checks[0]);
    state.writes = [];
    await reconcile(api, { apply: true, detailsUrl: 'https://example.org/runs/2' });
    assert.deepEqual(state.writes, []);
    assert.deepEqual(state.checks, [original]);
  }
});

test('dismissed and restored approvals replace completed checks', async () => {
  const { api, state } = server([pull(645, { labels: [label] })]);
  state.reviews.set(645, [review(1, specialist)]);
  await reconcile(api, { apply: true });
  state.reviews.set(645, [review(1, specialist, 'DISMISSED')]);
  await reconcile(api, { apply: true });
  assert.deepEqual(state.checks.map(check => check.conclusion), ['success', 'failure']);
  state.reviews.set(645, [review(1, specialist, 'DISMISSED'), review(2, specialist)]);
  await reconcile(api, { apply: true });
  assert.deepEqual(state.checks.map(check => check.conclusion), ['success', 'failure', 'success']);
});

test('changed explanations are published even when the conclusion stays the same', async () => {
  const { api, state } = server();
  state.reviews.set(645, []);
  await reconcile(api, { apply: true });
  state.reviews.set(645, [review(1, specialist)]);
  await reconcile(api, { apply: true });
  assert.deepEqual(state.checks.map(check => check.conclusion), ['failure', 'failure']);
  assert.match(state.checks[0].output.summary, /0\/2/);
  assert.match(state.checks[1].output.summary, /1\/2/);
});

test('unfinished checks can be completed without creating another run', async () => {
  for (const status of ['queued', 'in_progress']) {
    const { api, state } = server();
    state.checks.push({ id: 42, head_sha: head, name: 'review-policy', app: { id: 5288861 }, status });
    await reconcile(api, { apply: true });
    assert.equal(state.checks.length, 1);
    assert.equal(state.checks[0].conclusion, 'success');
    assert.equal(state.writes.every(write => write.method === 'PATCH' && write.path.endsWith('/42')), true);
  }
});

test('repeated evaluation errors preserve an identical completed failure', async () => {
  const { api, state } = server();
  state.before = ({ path }) => path === `${teamPath}/members` ? Response.json({}, { status: 403 }) : undefined;
  await assert.rejects(reconcile(api, { apply: true }), /403/);
  assert.equal(state.checks.at(-1).conclusion, 'failure');
  state.writes = [];
  await assert.rejects(reconcile(api, { apply: true }), /403/);
  assert.deepEqual(state.writes, []);
  assert.equal(state.checks.length, 1);
});

test('check lookup failures still publish failure for every known head', async () => {
  const secondHead = 'c'.repeat(40);
  const { api, state } = server([pull(), pull(646, { head: { sha: secondHead } })]);
  await reconcile(api, { apply: true });
  state.before = ({ path }) => path.includes('/commits/') ? Response.json({}, { status: 503 }) : undefined;
  await assert.rejects(reconcile(api, { apply: true }), /503/);
  for (const sha of [head, secondHead]) {
    assert.equal(state.checks.filter(check => check.head_sha === sha).at(-1).conclusion, 'failure');
  }
});

test('a publication error replaces successes already published for other heads', async () => {
  const secondHead = 'c'.repeat(40);
  const { api, state } = server([pull(), pull(646, { head: { sha: secondHead } })]);
  let failed = false;
  state.before = ({ path, body }) => {
    const sha = body?.head_sha ?? state.checks.find(check => path.endsWith(`/${check.id}`))?.head_sha;
    if (!failed && sha === secondHead && body?.conclusion === 'success') {
      failed = true;
      return Response.json({}, { status: 502 });
    }
  };
  await assert.rejects(reconcile(api, { apply: true }), /502/);
  assert.equal(failed, true);
  for (const sha of [head, secondHead]) {
    assert.equal(state.checks.filter(check => check.head_sha === sha).at(-1).conclusion, 'failure');
  }
});

test('composite action changes require two approvals, including deletions and renames', async () => {
  const files = [
    { filename: '.github/actions/test-ubuntu/action.yml', status: 'modified' },
    { filename: '.github/actions/test-macos/action.yml', status: 'removed' },
    { filename: '.github/actions/test-ubuntu/helper.sh', status: 'added' },
    { filename: 'scripts/helper.sh', previous_filename: '.github/actions/test-ubuntu/helper.sh', status: 'renamed' },
    { filename: '.github/actions/test-macos/helper.sh', previous_filename: 'scripts/helper.sh', status: 'renamed' },
  ];
  for (const file of files) {
    const { api, state } = server([pull(645, { labels: [label] })]);
    state.files.set(645, [file]);
    state.reviews.set(645, [review(1, specialist)]);
    const [result] = await reconcile(api);
    assert.equal(result.required, 2, JSON.stringify(file));
    assert.equal(result.passed, false);
    state.reviews.set(645, [review(1, specialist), review(2, maintainer)]);
    assert.equal((await reconcile(api))[0].passed, true);
  }
});

test('checks from a different app are never reused', async () => {
  const { api, state } = server();
  state.checks.push({ id: 42, head_sha: head, name: 'review-policy', app: { id: 15368 }, status: 'completed', conclusion: 'success' });
  await reconcile(api, { apply: true });
  assert.equal(state.checks.length, 2);
  assert.equal(state.checks[1].app.id, 5288861);
  assert.equal(state.writes.some(write => write.path.endsWith('/42')), false);
});

test('pagination includes later reviews, timeline events and team members', async () => {
  const { api, state } = server([pull(645, { labels: [label] })]);
  state.members = [...Array.from({ length: 100 }, (_, id) => ({ id: id + 1000 })), author, specialist];
  state.reviews.set(645, [...Array.from({ length: 100 }, (_, id) => review(id + 1, maintainer, 'COMMENTED')), review(101, specialist)]);
  state.timelines.set(645, [...Array.from({ length: 100 }, (_, id) => ({ id, event: 'commented' })), ...oop().timeline]);
  const results = await reconcile(api);
  assert.equal(results[0].passed, true);
  assert.equal(results[0].required, 1);
});

test('head or review changes during evaluation cannot publish stale success', async () => {
  for (const change of ['head', 'reviews']) {
    const { api, state } = server();
    let reads = 0;
    state.before = ({ path }) => {
      if (path === `${repoPath}/pulls/645/reviews` && ++reads === 2) {
        if (change === 'head') state.prs[0].head.sha = 'c'.repeat(40);
        else state.reviews.set(645, [review(1, specialist)]);
      }
    };
    await reconcile(api, { apply: true });
    assert.equal(state.checks[0].conclusion, 'failure', change);
  }
});

test('a newly opened PR sharing the commit invalidates a pending success', async () => {
  const { api, state } = server();
  let lists = 0;
  state.before = ({ path }) => {
    if (path === `${repoPath}/pulls` && ++lists === 2) state.prs.push(pull(646));
  };
  await reconcile(api, { apply: true });
  assert.equal(state.checks[0].conclusion, 'failure');
});

test('team lookup must resolve to the configured team identity', async () => {
  const { api, state } = server();
  state.teamId = 999;
  await assert.rejects(reconcile(api, { apply: true }), /team/i);
  assert.equal(state.checks[0].conclusion, 'failure');
});

test('API requests refuse external URLs and redirects', async () => {
  let calls = 0;
  const api = new GitHub('test-token', async (url, options) => {
    calls++;
    assert.equal(options.redirect, 'error');
    return Response.json({});
  });
  await assert.rejects(api.request('https://example.org/collect'), /path/i);
  assert.equal(calls, 0);
  await api.request('/orgs/SimVascular');
  assert.equal(calls, 1);
});
