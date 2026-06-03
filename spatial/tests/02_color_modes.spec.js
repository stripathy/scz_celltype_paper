import { test, expect } from '@playwright/test';

test('5 stable color modes switch cleanly without errors', async ({ page }) => {
  const errors = [];
  page.on('pageerror', e => errors.push(e.toString()));

  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  // SCZ has subclass / supertype / class / layer / depth as stable modes
  const stableModes = ['subclass', 'supertype', 'class', 'layer', 'depth'];
  for (const mode of stableModes) {
    await page.evaluate(m => {
      const btn = [...document.querySelectorAll('.mode-btn')].find(b => b.dataset.mode === m);
      btn.click();
    }, mode);
    await page.waitForFunction(m => window.colorMode === m, mode);
    await page.waitForTimeout(100);
  }

  expect(errors).toEqual([]);
});

test('confidence + margin modes appear when QC details toggled on', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  // QC-only modes should be hidden by default
  const initiallyHidden = await page.evaluate(() => {
    const conf = [...document.querySelectorAll('.mode-btn')].find(b => b.dataset.mode === 'confidence');
    return conf?.style.display === 'none';
  });
  expect(initiallyHidden).toBe(true);

  // Toggle Cell QC details on
  await page.evaluate(() => {
    const cb = document.getElementById('show-qc-details-toggle');
    cb.checked = true;
    cb.onchange({ target: cb });
  });
  await page.waitForFunction(() => window.showQcDetails === true);

  // Confidence + margin should now be visible
  const nowVisible = await page.evaluate(() => {
    const conf = [...document.querySelectorAll('.mode-btn')].find(b => b.dataset.mode === 'confidence');
    const margin = [...document.querySelectorAll('.mode-btn')].find(b => b.dataset.mode === 'margin');
    return conf?.style.display !== 'none' && margin?.style.display !== 'none';
  });
  expect(nowVisible).toBe(true);

  // Toggling off while in confidence should fall back to subclass
  await page.evaluate(() => {
    [...document.querySelectorAll('.mode-btn')].find(b => b.dataset.mode === 'confidence').click();
  });
  await page.waitForFunction(() => window.colorMode === 'confidence');
  await page.evaluate(() => {
    const cb = document.getElementById('show-qc-details-toggle');
    cb.checked = false;
    cb.onchange({ target: cb });
  });
  await page.waitForFunction(() => window.colorMode === 'subclass');
});
