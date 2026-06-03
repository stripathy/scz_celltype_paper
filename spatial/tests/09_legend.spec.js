import { test, expect } from '@playwright/test';

test('cell-type legend at ≤10 active, hides at >10 + in continuous modes', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  // All active → no per-type legend
  await page.evaluate(() => {
    [...document.querySelectorAll('#celltype-section .filter-bar button')]
      .find(b => b.textContent.trim() === 'All').click();
  });
  await page.waitForFunction(() => window.activeTypes.size > 10);
  const allLegend = await page.evaluate(() =>
    document.getElementById('legend-overlay').innerHTML);
  expect(allLegend).not.toContain('Active cell types');

  // Activate 5 → legend appears
  await page.evaluate(() => {
    [...document.querySelectorAll('#celltype-section .filter-bar button')]
      .find(b => b.textContent.trim() === 'None').click();
  });
  await page.waitForFunction(() => window.activeTypes.size === 0);
  await page.evaluate(() => {
    const rows = [...document.querySelectorAll('#celltype-filter .ct-row')].slice(0, 5);
    for (const r of rows) {
      const cb = r.querySelector('input[type="checkbox"]');
      cb.checked = true; cb.onchange();
    }
  });
  await page.waitForFunction(() => window.activeTypes.size === 5);
  await page.waitForTimeout(150);
  const five = await page.evaluate(() =>
    document.getElementById('legend-overlay').innerHTML);
  expect(five).toContain('Active cell types');

  // Switch to depth (continuous) — per-type legend hides
  await page.evaluate(() => {
    [...document.querySelectorAll('.mode-btn')].find(b => b.dataset.mode === 'depth').click();
  });
  await page.waitForFunction(() => window.colorMode === 'depth');
  await page.waitForTimeout(150);
  const depth = await page.evaluate(() =>
    document.getElementById('legend-overlay').innerHTML);
  expect(depth).toContain('Cortical Depth');
  expect(depth).not.toContain('Active cell types');
});
