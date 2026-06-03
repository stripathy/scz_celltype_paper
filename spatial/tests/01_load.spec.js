import { test, expect } from '@playwright/test';

test('page loads with no console errors', async ({ page }) => {
  const errors = [];
  page.on('console', msg => msg.type() === 'error' && errors.push(msg.text()));
  page.on('pageerror', e => errors.push(e.toString()));

  await page.goto('/');
  await page.waitForFunction(() => window.sampleData && window.indexData, null,
    { timeout: 15000 });

  expect(errors).toEqual([]);
  // Sanity: the cell count is right
  const nCells = await page.evaluate(() => window.sampleData.n_cells);
  expect(nCells).toBeGreaterThan(50000);
});
