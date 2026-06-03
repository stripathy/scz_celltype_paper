import { test, expect } from '@playwright/test';

test('cell-type color picker: sidebar swatch interactive + persists', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  const tag = await page.evaluate(() => {
    const sw = document.querySelector('#celltype-filter .ct-row .ct-swatch');
    return { tag: sw?.tagName, type: sw?.type };
  });
  expect(tag.tag).toBe('INPUT');
  expect(tag.type).toBe('color');

  const targetName = await page.evaluate(() => {
    const row = document.querySelector('#celltype-filter .ct-row');
    const sw = row.querySelector('.ct-swatch');
    const lbl = row.querySelector('.ct-label').textContent;
    sw.value = '#ff00ff';
    sw.dispatchEvent(new Event('input', { bubbles: true }));
    return lbl;
  });

  await page.waitForFunction(name =>
    window.customCellTypeColors[window.colorMode]?.[name] === '#ff00ff',
    targetName);
});

test('cell-type color override is mode-scoped (subclass vs supertype)', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  // Set a subclass color override
  const subName = await page.evaluate(() => {
    const row = document.querySelector('#celltype-filter .ct-row');
    const sw = row.querySelector('.ct-swatch');
    sw.value = '#abcdef';
    sw.dispatchEvent(new Event('input', { bubbles: true }));
    return row.querySelector('.ct-label').textContent;
  });
  await page.waitForFunction(n => window.customCellTypeColors.subclass?.[n] === '#abcdef', subName);

  // Switch to supertype mode — subclass override should NOT bleed
  await page.evaluate(() => {
    [...document.querySelectorAll('.mode-btn')].find(b => b.dataset.mode === 'supertype').click();
  });
  await page.waitForFunction(() => window.colorMode === 'supertype');

  // No supertype override yet
  const supScopedEmpty = await page.evaluate(() =>
    Object.keys(window.customCellTypeColors.supertype || {}).length === 0);
  expect(supScopedEmpty).toBe(true);

  // Switch back to subclass — original override still there
  await page.evaluate(() => {
    [...document.querySelectorAll('.mode-btn')].find(b => b.dataset.mode === 'subclass').click();
  });
  await page.waitForFunction(() => window.colorMode === 'subclass');
  const stillThere = await page.evaluate(n =>
    window.customCellTypeColors.subclass?.[n] === '#abcdef', subName);
  expect(stillThere).toBe(true);
});
