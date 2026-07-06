"""Capture headless screenshots of the live POLARIS site (served locally) into docs/screenshots/."""
import pathlib, sys
from playwright.sync_api import sync_playwright

OUT = pathlib.Path("docs/screenshots"); OUT.mkdir(parents=True, exist_ok=True)
URL = "http://localhost:8753/"


def run():
    with sync_playwright() as p:
        b = p.chromium.launch()
        # ---- desktop ----
        pg = b.new_page(viewport={"width": 1400, "height": 900}, device_scale_factor=1.5)
        pg.goto(URL, wait_until="networkidle"); pg.wait_for_timeout(1600)
        pg.screenshot(path=str(OUT / "01_hero.png"))

        def shot(sel, name, settle=750):
            el = pg.query_selector(sel)
            el.scroll_into_view_if_needed(); pg.wait_for_timeout(settle)
            pg.eval_on_selector("#nav", "n=>n.style.visibility='hidden'")  # avoid sticky-nav bleed
            el.screenshot(path=str(OUT / name))
            pg.eval_on_selector("#nav", "n=>n.style.visibility='visible'")

        shot("#guardrails", "02_principles.png")
        shot("#pipeline", "03_pipeline.png")
        # explorer with the FTO disagreement selected
        pg.eval_on_selector("#explorer", "e=>e.scrollIntoView()"); pg.wait_for_timeout(500)
        try:
            pg.click('.exp-row[data-loc="FTO"]', timeout=4000); pg.wait_for_timeout(900)
        except Exception as e:
            print("FTO click skipped:", e)
        shot("#explorer", "04_explorer.png")
        shot("#vignette", "05_fto_vignette.png", 1400)
        shot("#validation", "06_validation.png", 1000)
        shot("#generalization", "07_generalization.png")
        shot("#honesty", "08_honesty.png")

        # ---- mobile ----
        mp = b.new_page(viewport={"width": 390, "height": 844}, device_scale_factor=2)
        mp.goto(URL, wait_until="networkidle"); mp.wait_for_timeout(1500)
        mp.screenshot(path=str(OUT / "09_mobile_hero.png"))
        mp.eval_on_selector("#explorer", "e=>e.scrollIntoView()"); mp.wait_for_timeout(1200)
        mp.screenshot(path=str(OUT / "10_mobile_explorer.png"))
        b.close()
    print("saved:", sorted(p.name for p in OUT.glob("*.png")))


if __name__ == "__main__":
    run()
