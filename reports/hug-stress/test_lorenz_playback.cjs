// Run with Node and Playwright installed: node test_lorenz_playback.cjs [fragment.html]
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const {pathToFileURL} = require('node:url');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const target = path.resolve(process.argv[2] || path.join(__dirname, 'lorenz-hug-fragment.html'));
const records = JSON.parse(fs.readFileSync(path.join(__dirname, 'lorenz-data/presets-matched.json'), 'utf8'));
let browser;
(async () => {
  browser = await chromium.launch({channel:'chrome', headless:true});
  const errors = [], checks = [];
  async function open(state) {
    const page = await browser.newPage();
    page.on('pageerror', error => errors.push(String(error)));
    await page.addInitScript(initial => {
      window.__saves = [];
      window.openai = {widgetState:initial, setWidgetState(value) {
        window.__saves.push(value); return Promise.resolve();
      }};
      // Deliberately separate the animation clock from performance.now().
      Object.defineProperty(performance, 'now', {value:() => 1000000});
      let next = 0;
      const frames = new Map();
      window.requestAnimationFrame = callback => {frames.set(++next, callback); return next;};
      window.cancelAnimationFrame = id => frames.delete(id);
      window.__step = time => {
        const pending = [...frames.values()]; frames.clear();
        for (const callback of pending) callback(time);
      };
    }, state);
    await page.goto(pathToFileURL(target).href);
    return page;
  }
  const snapshot = page => page.evaluate(() => document.getElementById('lorenz-hug-lab').lorenzSnapshot());
  async function check(page, key, index, running = false) {
    const state = await snapshot(page);
    assert.equal(state.key, key); assert.equal(state.index, index); assert.equal(state.running, running);
    assert.deepEqual(state.sample, records.find(r => r.key === key).samples[index]);
  }
  const page = await open(null);
  await check(page, 'chaotic', 0);
  assert.equal(await page.evaluate(() => __saves.length), 0);
  await page.locator('#lh-play').click();
  await page.evaluate(() => __step(0)); // Previously throws: sample is not iterable.
  await check(page, 'chaotic', 0, true);
  await page.evaluate(() => __step(80));
  await check(page, 'chaotic', 2, true);
  await page.evaluate(() => __step(-10));
  await check(page, 'chaotic', 2, true);
  await page.locator('#lh-play').click();
  await page.evaluate(() => __step(2000));
  await check(page, 'chaotic', 2);
  await page.locator('#lh-play').click();
  await page.evaluate(() => {__step(100); __step(180);});
  await check(page, 'chaotic', 4, true);
  await page.evaluate(() => __step(100000));
  await check(page, 'chaotic', 400);
  await page.locator('#lh-play').click();
  await page.evaluate(() => __step(100));
  await check(page, 'chaotic', 0, true);
  checks.push('clock mismatch, backward timestamp, pause, resume, end and replay');
  // A delayed state event must stop playback without saving it back to the host.
  const before = await page.evaluate(() => __saves.length);
  await page.evaluate(() => window.dispatchEvent(new CustomEvent('openai:set_globals', {
    detail:{globals:{widgetState:{privateContent:{key:'opening', index:123}}}}
  })));
  await page.evaluate(() => __step(1000));
  await check(page, 'opening', 123);
  assert.equal(await page.evaluate(() => __saves.length), before);
  checks.push('delayed saved state cancels playback without a save loop');
  for (const record of records) {
    await page.locator('#lh-mode').selectOption(record.key);
    for (const index of [0, 10, 200, 400]) {
      await page.locator('#lh-time').evaluate((el, i) => {
        el.value = i; el.dispatchEvent(new Event('input')); el.dispatchEvent(new Event('change'));
      }, index);
      await check(page, record.key, index);
    }
  }
  checks.push('all five presets and 20 slider selections match recorded samples');
  await page.emulateMedia({reducedMotion:'reduce'});
  await page.locator('#lh-play').click();
  await check(page, 'opening', 400);
  checks.push('reduced motion');
  await page.close();
  for (const [value, expected] of [
    [undefined, 0], [null, 0], ['invalid', 0], [{}, 0], [[], 0], [true, 0],
    [NaN, 0], [Infinity, 0], [-Infinity, 0], [-10, 0], [999, 400], [12.6, 13], [200, 200]
  ]) {
    const restored = await open({privateContent:{key:'steady', index:value}});
    await check(restored, 'steady', expected);
    assert.equal(await restored.evaluate(() => __saves.length), 0);
    await restored.close();
  }
  const incompatible = await open({privateContent:{key:'another-visual', index:150}});
  await check(incompatible, 'chaotic', 0);
  await incompatible.close();
  checks.push('13 saved index cases and incompatible preset fallback');
  assert.deepEqual(errors, []);
  const result = {checks, errors};
  fs.writeFileSync(path.join(__dirname, 'lorenz-data/playback-regression.json'), JSON.stringify(result, null, 2) + '\n');
  console.log(JSON.stringify(result));
  await browser.close();
})().catch(async error => {console.error(error); if (browser) await browser.close(); process.exitCode = 1;});
