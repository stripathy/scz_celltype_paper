import { test, expect } from '@playwright/test';

test('cell-type search filters rows', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  // Type "Pvalb" into the search box
  await page.fill('#celltype-search', 'Pvalb');
  // After input, the filter rebuilds — wait for the rows to settle
  await page.waitForFunction(() => {
    const labels = [...document.querySelectorAll('#celltype-filter .ct-row .ct-label')]
      .map(el => el.textContent);
    return labels.length > 0 && labels.every(l => l.toLowerCase().includes('pvalb'));
  });

  const labels = await page.$$eval('#celltype-filter .ct-row .ct-label',
    els => els.map(el => el.textContent));
  expect(labels.length).toBeGreaterThan(0);
  expect(labels.every(l => l.toLowerCase().includes('pvalb'))).toBe(true);

  // Clear search — should restore all rows (SCZ default subclass ~25 rows)
  await page.fill('#celltype-search', '');
  await page.waitForFunction(() =>
    document.querySelectorAll('#celltype-filter .ct-row').length > 5);
});
