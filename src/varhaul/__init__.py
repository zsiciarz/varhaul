"""List eclipsing binaries with a primary or secondary minimum tonight.

Reads linear elements in the Krakow ephemeris format
(``CON VAR  ALL|PRI|SEC  M0 (err)  P (err) ...``, M0 in HJD) and takes each
star's position, magnitude range and type from the GCVS catalogue
(``gcvs5.txt``). Writes the minima falling in the night's dark window while
the star is above ``--min-alt`` to a static HTML page, as a table sorted by
time. Minima above ``--max-alt`` (too close to the zenith for an alt-az
mount) are kept but flagged. Secondary minima come from a SEC line when the star has one, else
from phase 0.5 of the primary elements.

Needs only astropy:

    uv run varhaul [DATE] [--lat 52.23 --lon 21.01 --min-alt 20] [-o FILE]
"""

import argparse
import datetime as dt
import html
import math
import re
import sys
from importlib.resources import files
from pathlib import Path
from string import Template
from zoneinfo import ZoneInfo

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[2]  # repo root (src/varhaul/__init__.py)
WARSAW = (52.2297, 21.0122)
ROW = re.compile(
    r"^(?P<con>[A-Z]{3}) (?P<var>\S+(?: \d)?)\s+(?P<kind>ALL|PRI|SEC|all|pri|sec)\s+"
    r"(?P<m0>\d+\.\d+)\s*\((?P<m0err>[\d*]+)\)\s+(?P<p>\d+\.\d+)\s*\((?P<perr>[\d*]+)\)"
)
# GCVS J2000 position, e.g. ``231110.10  +530133.0`` (seconds may end in ``.``).
COORDS = re.compile(
    r"(\d\d)(\d\d)(\d\d(?:\.\d*)?)\s*([+-])(\d\d)(\d\d)(\d\d(?:\.\d*)?)"
)
# Ephemeris names that GCVS lists under another designation.
ALIASES = {"BOO 44": "BOO I"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "date",
        nargs="?",
        type=dt.date.fromisoformat,
        default=dt.date.today(),
        help="evening date of the night (default: today)",
    )
    parser.add_argument(
        "--lat", type=float, default=WARSAW[0], help="latitude, degrees N"
    )
    parser.add_argument(
        "--lon", type=float, default=WARSAW[1], help="longitude, degrees E"
    )
    parser.add_argument(
        "--tz",
        type=ZoneInfo,
        default=ZoneInfo("Europe/Warsaw"),
        help="timezone for times",
    )
    parser.add_argument(
        "--min-alt",
        type=float,
        default=20.0,
        help="horizon altitude, degrees (default: 20)",
    )
    parser.add_argument(
        "--max-alt",
        type=float,
        default=64.0,
        help="flag minima above this altitude, degrees (default: 64)",
    )
    parser.add_argument(
        "--sun-alt",
        type=float,
        default=-12.0,
        help="Sun altitude ending twilight (default: -12)",
    )
    parser.add_argument(
        "--ephem", type=Path, default=ROOT / "data/ephem.txt", help="ephemeris file"
    )
    parser.add_argument(
        "--gcvs", type=Path, default=ROOT / "data/gcvs5.txt", help="GCVS catalogue"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="HTML file to write (default: minima-DATE.html)",
    )
    args = parser.parse_args()
    if args.output is None:
        args.output = Path(f"minima-{args.date}.html")
    for path in (args.ephem, args.gcvs):
        if not path.is_file():
            parser.error(f"file not found: {path}")
    return args


def err(value: str, digits: str) -> float | None:
    """An error given in units of the value's last digit, or None if unknown."""
    if digits == "*":
        return None
    decimals = len(value.split(".")[1])
    return int(digits) * 10**-decimals


def key(con: str, var: str) -> str:
    """A normalised star name shared by both catalogues: GCVS ``V0683 Lac``
    and ephemeris ``LAC V683`` -> ``LAC V683``; ``mu. 1 Sco`` -> ``SCO MU1``."""
    var = re.sub(r"[\s.]", "", var).upper()
    if m := re.fullmatch(r"V0*(\d+)", var):
        var = f"V{m.group(1)}"
    return f"{con.upper()} {var}"


def read_ephemeris(path: Path) -> dict[str, dict[str, dict]]:
    ephem: dict[str, dict[str, dict]] = {}
    for line in path.read_text().splitlines():
        if m := ROW.match(line):
            ephem.setdefault(key(m["con"], m["var"]), {})[m["kind"].upper()] = {
                "m0": float(m["m0"]),
                "m0err": err(m["m0"], m["m0err"]),
                "p": float(m["p"]),
                "perr": err(m["p"], m["perr"]),
            }
    return ephem


def read_gcvs(path: Path) -> dict[str, dict]:
    """GCVS stars with a known position, by ``key()``."""
    stars = {}
    for line in path.read_text(encoding="latin-1").splitlines():
        f = line.split("|")
        if len(f) < 8 or not (c := COORDS.match(f[2].strip())):
            continue
        h, m, s, sign, d, dm, ds = c.groups()
        name = " ".join(f[1].replace("*", "").split())
        var, con = name.rsplit(" ", 1)
        stars[key(con, var)] = {
            "name": name,
            "ra": int(h) + int(m) / 60 + float(s) / 3600,
            "dec": (1 if sign == "+" else -1)
            * (int(d) + int(dm) / 60 + float(ds) / 3600),
            "type": f[3].strip(),
            "max_mag": re.sub(r"\s", "", f[4]),
            "min_mag": re.sub(r"\s", "", f[5]),
            "band": f[7].strip(),
        }
    return stars


def lookup(gcvs: dict[str, dict], k: str) -> dict | None:
    """The GCVS entry for an ephemeris key, allowing for aliases and
    component suffixes (``HER V994A`` -> ``HER V994``, shown as ``V0994 Her A``)."""
    if star := gcvs.get(ALIASES.get(k, k)):
        return star
    if (m := re.fullmatch(r"(.* V\d+)([A-Z])", k)) and (star := gcvs.get(m[1])):
        return {**star, "name": f"{star['name']} {m[2]}"}
    return None


def night_window(
    date: dt.date, location: EarthLocation, tz: ZoneInfo, sun_alt: float
) -> tuple[Time, Time]:
    """The longest stretch the Sun stays below ``sun_alt`` in the 24 hours
    from local noon on ``date``, sampled every 2 minutes."""
    noon = dt.datetime.combine(date, dt.time(12), tzinfo=tz).astimezone(dt.UTC)
    times = Time(noon) + np.arange(0, 24 * 60, 2) * u.min
    alt = get_sun(times).transform_to(AltAz(obstime=times, location=location)).alt.deg
    dark = np.concatenate([[False], alt < sun_alt, [False]])
    edges = np.flatnonzero(np.diff(dark.astype(int)))
    if not len(edges):
        sys.exit(f"The Sun never gets below {sun_alt:g}° on {date}.")
    starts, ends = edges[::2], edges[1::2] - 1
    longest = np.argmax(ends - starts)
    return times[starts[longest]], times[ends[longest]]


def minima(el: dict[str, dict], ltt: float, start: Time, end: Time) -> list[dict]:
    """A star's minima between ``start`` and ``end`` (no altitude check).
    ``ltt`` is the heliocentric light-time correction in days."""
    pri = el.get("PRI") or el.get("ALL")
    if pri is None:  # only a secondary ephemeris
        return []
    sec = el.get("SEC") or {**pri, "m0": pri["m0"] + pri["p"] / 2}
    # HJD -> UTC JD: subtract the heliocentric light-time correction.
    lo, hi = start.jd + ltt, end.jd + ltt
    found = []
    for kind, e in (("I", pri), ("II", sec)):
        n = math.ceil((lo - e["m0"]) / e["p"])
        while (hjd := e["m0"] + n * e["p"]) <= hi:
            sigma = None
            if e["m0err"] is not None and e["perr"] is not None:
                sigma = math.hypot(e["m0err"], n * e["perr"]) * 24 * 60
            found.append({"jd": hjd - ltt, "min": kind, "sigma": sigma})
            n += 1
    return found


COMPASS = "N NNE NE ENE E ESE SE SSE S SSW SW WSW W WNW NW NNW".split()
MAG = re.compile(r"-?\d+(?:\.\d*)?")


def mag_range(star: dict) -> tuple[float | None, float | None]:
    """A star's (brightest, faintest) magnitude from its GCVS fields.
    ``:`` (uncertain) and ``<``/``>`` (limits) are ignored, so ``<16.`` counts
    as 16; a bracketed minimum like ``(0.29)`` is an amplitude."""
    bright = float(m[0]) if (m := MAG.search(star["max_mag"])) else None
    faint = float(m[0]) if (m := MAG.search(star["min_mag"])) else None
    if star["min_mag"].startswith("("):
        faint = None if bright is None or faint is None else bright + faint
    return bright, faint


def render(
    args: argparse.Namespace,
    start: Time,
    end: Time,
    rows: list[dict],
    unmatched: list[str],
    local,
) -> str:
    """The night's minima as a self-contained HTML page."""
    cells = []
    for r in rows:
        star, el = r["star"], r["el"]
        period = (el.get("PRI") or el["ALL"])["p"] * 24
        sigma = "?" if r["sigma"] is None else f"{r['sigma']:.0f}"
        sigma_sort = "" if r["sigma"] is None else f"{r['sigma']:.1f}"
        compass = COMPASS[round(r["az"] / 22.5) % 16]
        cls = ' class="zenith"' if r["alt"] > args.max_alt else ""
        bright, faint = mag_range(star)
        mags = "".join(
            f' data-{k}="{v:g}"' for k, v in (("bright", bright), ("faint", faint))
            if v is not None
        )
        e = html.escape
        cells.append(
            f'<tr{cls}{mags}><td data-sort="{r["jd"]:.5f}">{r["time"]:%H:%M}</td>'
            f'<td class="star">{e(star["name"])}</td>'
            f"<td>{r['min']}</td><td>{e(star['type'])}</td>"
            f'<td data-sort="{"" if bright is None else f"{bright:g}"}">'
            f"{e(star['max_mag'])}–{e(star['min_mag'])} {e(star['band'])}</td>"
            f'<td class="num">{period:.1f}</td>'
            f'<td class="num alt">{r["alt"]:.0f}°</td>'
            f'<td class="num" data-sort="{r["az"]:.1f}">{r["az"]:.0f}° {compass}</td>'
            f'<td class="num" data-sort="{sigma_sort}">{sigma}</td></tr>'
        )
    template = Template(files(__package__).joinpath("template.html").read_text("utf-8"))
    return template.substitute(
        date=args.date.isoformat(),
        location=f"{args.lat:.4f}°, {args.lon:.4f}° ({html.escape(str(args.tz))})",
        sun_alt=f"{args.sun_alt:g}",
        window=f"{local(start):%Y-%m-%d %H:%M} – {local(end):%H:%M}",
        min_alt=f"{args.min_alt:g}",
        max_alt=f"{args.max_alt:g}",
        count=len(rows),
        rows="\n".join(cells),
        unmatched_count=len(unmatched),
        unmatched=html.escape(", ".join(unmatched)) or "none",
        generated=dt.datetime.now(args.tz).strftime("%Y-%m-%d %H:%M %Z"),
    )


def main() -> None:
    args = parse_args()
    location = EarthLocation(lat=args.lat * u.deg, lon=args.lon * u.deg)
    start, end = night_window(args.date, location, args.tz, args.sun_alt)
    local = lambda t: t.to_datetime(dt.UTC).astimezone(args.tz)  # noqa: E731
    ephem = read_ephemeris(args.ephem)
    gcvs = read_gcvs(args.gcvs)
    stars, unmatched = [], []
    for k, el in ephem.items():
        if star := lookup(gcvs, k):
            stars.append((star, el))
        else:
            unmatched.append(k)

    coords = SkyCoord(
        [s["ra"] for s, _ in stars] * u.hourangle, [s["dec"] for s, _ in stars] * u.deg
    )
    ltts = (
        (start + (end - start) / 2)
        .light_travel_time(coords, "heliocentric", location)
        .jd
    )
    found = []
    for i, ((star, el), ltt) in enumerate(zip(stars, ltts)):
        found += [{**m, "i": i, "star": star, "el": el} for m in minima(el, ltt, start, end)]

    rows = []
    if found:
        times = Time([m["jd"] for m in found], format="jd", scale="utc")
        altaz = coords[[m["i"] for m in found]].transform_to(
            AltAz(obstime=times, location=location)
        )
        for m, t, alt, az in zip(found, times, altaz.alt.deg, altaz.az.deg):
            if alt >= args.min_alt:
                rows.append({**m, "time": local(t), "alt": alt, "az": az})
    rows.sort(key=lambda r: r["time"])

    page = render(args, start, end, rows, unmatched, local)
    args.output.write_text(page, encoding="utf-8")
    print(
        f"{args.output}: {len(rows)} minima, "
        f"{local(start):%H:%M}–{local(end):%H:%M} dark window; "
        f"not in GCVS: {', '.join(unmatched) or 'none'}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
