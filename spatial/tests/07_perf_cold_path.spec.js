import { test, expect } from '@playwright/test';

// Regression test for the 4.6s freeze caused by Path2D cold-path tessellation.
// CRITICAL: do NOT call selectAllTypes() to warm up — the warm path was always
// fast; the freeze only appears on the very first state-changing click.
test('first None click after page load is under 200ms', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);
  // Make sure render has happened once with the default state
  await page.waitForFunction(() => window.viewScale > 0);

  const ms = await page.evaluate(async () => {
    const noneBtn = [...document.querySelectorAll('#celltype-section .filter-bar button')]
      .find(b => b.textContent.trim() === 'None');
    const t0 = performance.now();
    noneBtn.click();
    // Force a paint to settle
    await new Promise(r => requestAnimationFrame(r));
    return performance.now() - t0;
  });

  // Original bug: 4662ms. After fix: 12ms. Set ceiling at 200ms for headroom.
  expect(ms).toBeLessThan(200);
});

test('render budget under 50ms at low zoom with show-deselected on', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  // Filter to one type to maximize the dim-layer load
  await page.click('#solo-btn');
  await page.evaluate(() => {
    document.querySelector('#celltype-filter .ct-row').click();
  });
  await page.waitForFunction(() => window.activeTypes.size === 1);
  await page.waitForFunction(() => window.showDeselectedCells === true);

  const ms = await page.evaluate(() => {
    // Median of 5 renders
    const ts = [];
    for (let i = 0; i < 5; i++) {
      const t0 = performance.now();
      window.render();
      ts.push(performance.now() - t0);
    }
    ts.sort((a, b) => a - b);
    return ts[2];
  });

  expect(ms).toBeLessThan(50);
});
