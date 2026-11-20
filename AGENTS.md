# AGENTS.md

varhaul plans eclipsing-binary observations: for a date and location it writes a static HTML page of all primary/secondary minima in that night's dark window. User-facing docs are in README.md.

## Commands

```sh
uv sync                                  # Python 3.14, only dependency is astropy
uv run varhaul [DATE] -o /tmp/out.html   # writes the page; one-line summary on stderr
```

There are no tests or linter. To check a change, run the tool and compare the stderr summary (minima count, unmatched stars) with a run from before the change.

## Layout

- `src/varhaul/__init__.py`: the whole pipeline, driven by `main()`.
- `src/varhaul/template.html`: the page, filled in by `render()` via `string.Template`. Write a literal `$` as `$$`. The page must stay self-contained (inline CSS/JS, no external requests).
- `data/ephem.txt` (Suhora elements) and `data/gcvs5.txt` (GCVS positions, magnitudes, types): third-party catalogues, git-ignored. **Never commit or redistribute them.** Default paths resolve from `__file__`, which relies on the editable install `uv sync` does.

## Gotchas

- Ephemeris epochs are HJD; convert with the heliocentric light-travel time before comparing to UTC.
- Both catalogues name stars differently. Always match through `key()` (`V0683 Lac` and `LAC V683` → `LAC V683`). Special cases go in `ALIASES`. `CEN PSI` and `HER 68` have no GCVS entry, so they're expected in the unmatched list.
- Keep the astropy work vectorised over all ~3,800 stars (one `SkyCoord`, one `AltAz` transform); a full run takes about 2 s this way, and per-star astropy calls are far slower.
- Catalogue text goes into HTML, so pass it through `html.escape`.
