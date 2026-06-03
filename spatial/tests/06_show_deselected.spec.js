import { test, expect } from '@playwright/test';

test('Show deselected: defaults ON + dim count appears', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  const defaults = await page.evaluate(() => ({
    state: window.showDeselectedCells,
    checked: document.getElementById('show-deselected-toggle')?.checked,
  }));
  expect(defaults.state).toBe(true);
  expect(defaults.checked).toBe(true);

  // Filter to one type via Solo
  await page.click('#solo-btn');
  await page.evaluate(() => {
    document.querySelector('#celltype-filter .ct-row').click();
  });
  await page.waitForFunction(() => window.activeTypes.size === 1);

  const status = await page.evaluate(() => document.getElementById('status').textContent);
  expect(status).toMatch(/\(\+\s*[\d,]+\s*dimmed\)/);
});

test('hover over deselected cell shows [deselected] tooltip badge', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  await page.click('#solo-btn');
  await page.evaluate(() => {
    document.querySelector('#celltype-filter .ct-row').click();
  });
  await page.waitForFunction(() => window.activeTypes.size === 1);

  const dispatched = await page.evaluate(() => {
    // SCZ uses individual subclass/supertype/class/layer arrays, not a clusterings dict
    let cats, idx;
    const m = window.colorMode;
    if (m === 'supertype') { cats = window.sampleData.supertype_cats; idx = window.sampleData.supertype; }
    else if (m === 'class') { cats = window.sampleData.class_cats; idx = window.sampleData.class; }
    else if (m === 'layer') { cats = window.sampleData.layer_cats; idx = window.sampleData.layer; }
    else { cats = window.sampleData.subclass_cats; idx = window.sampleData.subclass; }
    let target = -1;
    for (let i = 0; i < window.sampleData.n_cells; i++) {
      // Skip qc-failed cells if hideQcFail is on
      if (window.hideQcFail && window.sampleData.qc_status?.[i] > 0) continue;
      const name = cats[idx[i]];
      if (!window.activeTypes.has(name)) { target = i; break; }
    }
    if (target < 0) return false;
    const sx = window.sampleData.x[target] * window.viewScale + window.viewX;
    const sy = window.sampleData.y[target] * window.viewScale + window.viewY;
    const canvas = document.querySelector('canvas');
    const rect = canvas.getBoundingClientRect();
    canvas.dispatchEvent(new MouseEvent('mousemove', {
      clientX: rect.left + sx, clientY: rect.top + sy, bubbles: true,
    }));
    return true;
  });
  expect(dispatched).toBe(true);

  await page.waitForFunction(() =>
    document.getElementById('tooltip')?.style.display === 'block',
    null, { timeout: 2000 });
  const tt = await page.evaluate(() => document.getElementById('tooltip').textContent);
  expect(tt).toContain('deselected');
});

test('hideQcFail interaction: dim count excludes qc-fail cells', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => window.sampleData);

  // Solo on first type to get a small active set
  await page.click('#solo-btn');
  await page.evaluate(() => {
    document.querySelector('#celltype-filter .ct-row').click();
  });
  await page.waitForFunction(() => window.activeTypes.size === 1);

  // Capture dim count with hideQcFail=true (default)
  const dimWithHide = await page.evaluate(() => {
    const s = document.getElementById('status').textContent;
    const m = s.match(/\(\+\s*([\d,]+)\s*dimmed\)/);
    return m ? parseInt(m[1].replace(/,/g, ''), 10) : null;
  });

  // Toggle hideQcFail off
  await page.evaluate(() => {
    const cb = document.getElementById('hide-qc-fail-toggle');
    cb.checked = false;
    cb.onchange({ target: cb });
  });
  await page.waitForFunction(() => window.hideQcFail === false);

  const dimWithoutHide = await page.evaluate(() => {
    const s = document.getElementById('status').textContent;
    const m = s.match(/\(\+\s*([\d,]+)\s*dimmed\)/);
    return m ? parseInt(m[1].replace(/,/g, ''), 10) : null;
  });

  // dimWithoutHide should be exactly the number of qc-fail cells more
  const qcFailCount = await page.evaluate(() => {
    let c = 0;
    for (let i = 0; i < window.sampleData.qc_status.length; i++) {
      if (window.sampleData.qc_status[i] > 0) c++;
    }
    return c;
  });

  expect(dimWithoutHide - dimWithHide).toBe(qcFailCount);
});
