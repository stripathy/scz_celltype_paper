import { test, expect } from '@playwright/test';

// Pixel-checksum regression test. The first run records a reference image;
// subsequent runs compare against it. Rebaseline with `npm run test:update`.
test('default canvas render matches reference', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);
  // Force a render with the default state and a fixed view
  await page.evaluate(() => {
    if (typeof window.fitView === 'function') window.fitView();
    if (typeof window.render === 'function') window.render();
  });
  await page.waitForTimeout(200);
  await expect(page.locator('canvas')).toHaveScreenshot('default-render.png',
    { maxDiffPixelRatio: 0.05 });
});
