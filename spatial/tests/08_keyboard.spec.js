import { test, expect } from '@playwright/test';

test('keyboard shortcuts: a / x / f', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  // 'a' selects all
  await page.evaluate(() => {
    [...document.querySelectorAll('#celltype-section .filter-bar button')]
      .find(b => b.textContent.trim() === 'None').click();
  });
  await page.waitForFunction(() => window.activeTypes.size === 0);
  await page.keyboard.press('a');
  await page.waitForFunction(() => window.activeTypes.size > 5);

  // 'x' toggles showDeselectedCells
  const before = await page.evaluate(() => window.showDeselectedCells);
  await page.keyboard.press('x');
  await page.waitForFunction(prev => window.showDeselectedCells !== prev, before);

  // 'f' toggles hideQcFail
  const beforeF = await page.evaluate(() => window.hideQcFail);
  await page.keyboard.press('f');
  await page.waitForFunction(prev => window.hideQcFail !== prev, beforeF);

  // Search input should suppress the shortcut
  await page.click('#celltype-search');
  await page.evaluate(() => {
    [...document.querySelectorAll('#celltype-section .filter-bar button')]
      .find(b => b.textContent.trim() === 'None').click();
  });
  await page.waitForFunction(() => window.activeTypes.size === 0);
  await page.keyboard.type('a');
  expect(await page.inputValue('#celltype-search')).toBe('a');
  expect(await page.evaluate(() => window.activeTypes.size)).toBe(0);
});
