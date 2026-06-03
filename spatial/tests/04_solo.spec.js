import { test, expect } from '@playwright/test';

test('Solo button + row click sets sole-active type', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  // Initial state: all types active
  const initialActive = await page.evaluate(() => window.activeTypes.size);
  // SCZ default subclass has ~25 types; supertype has ~140
  expect(initialActive).toBeGreaterThan(5);

  // Click Solo
  await page.click('#solo-btn');
  await page.waitForFunction(() => window.soloMode === true);

  // Click first row — should make it the sole active type
  const firstRowLabel = await page.evaluate(() => {
    const row = document.querySelector('#celltype-filter .ct-row');
    const lbl = row?.querySelector('.ct-label')?.textContent;
    row.click();
    return lbl;
  });
  await page.waitForFunction(() => window.activeTypes.size === 1);

  const soloType = await page.evaluate(() => window.soloType);
  expect(soloType).toBe(firstRowLabel);

  // Click Solo again to exit
  await page.click('#solo-btn');
  await page.waitForFunction(() => window.soloMode === false);
  // After exit, all types restored
  const finalActive = await page.evaluate(() => window.activeTypes.size);
  expect(finalActive).toBe(initialActive);
});
