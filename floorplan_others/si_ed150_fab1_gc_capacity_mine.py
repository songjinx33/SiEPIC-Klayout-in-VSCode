# %%
"""Real numbered Universal_GC loopbacks in a fixed ZEP layout area.

Each corner reserves one normal two-GC loopback and one three-GC Y-branch
alignment device. Its center GC connects to the splitter input; the upper
and lower GCs connect to the two outputs. All remaining legal positions
contain numbered normal loopbacks. Packing and final checks use real cell
bounding boxes, including cladding, routes and recognition layers.
Dimensions are micrometers; packing uses integer database units.
Capacity is for an axis-aligned regular grid.
"""
from dataclasses import dataclass
from pathlib import Path
import math

area_width = 12000
area_height = 6000
margin_x = 100
margin_y = 100
device_gap_x = 20
device_gap_y = 20
gc_pitch = 127.0
# Optional minimum footprint; None means measured component dimensions.
required_x_size = None
required_y_size = None
PDK_ROOT = Path(r"D:\Github\SiEPICfab-EBeam-ZEP-PDK\SiEPICfab_EBeam_ZEP")
tech_name = "SiEPICfab_EBeam_ZEP"  # actual installed technology spelling
GC_SETTINGS = dict(
    gc_cell_name="Universal_GC", gc_library=tech_name + "_Beta", gc_pitch=gc_pitch,
    parameters=dict(target_length=25.0, duty_cycle=0.45, pitch=0.7,
                    radius=25.0, y_span=15.0, L_extra=5.0,
                    wg_width=0.5, CladWidth=2.0),
)  # Geometry inputs; optical performance requires simulation/measurement.
WG_SETTINGS = {"wg_type": "Strip TE 1550 nm, w=500 nm (core-clad)"}
wg_bend_radius = 20
show_layout = True
BASE_DIR = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()


@dataclass(frozen=True)
class Rect:
    left: int
    bottom: int
    right: int
    top: int

    @property
    def width(self):
        return self.right - self.left

    @property
    def height(self):
        return self.top - self.bottom


def separated(a, b, gx, gy):
    """At least the requested gap on one axis; exact equality is allowed."""
    return (a.right + gx <= b.left or b.right + gx <= a.left
            or a.top + gy <= b.bottom or b.top + gy <= a.bottom)


def packing(width, height, mx, my, gx, gy, footprint, normal, ybranch):
    """Pure integer geometry shared by both codes, including corner exclusion."""
    values = (width, height, mx, my, gx, gy, *footprint, *normal, *ybranch)
    if any(v < 0 for v in values) or min(width, height, *footprint, *normal, *ybranch) <= 0:
        raise ValueError("Area and footprints must be positive; margins/gaps nonnegative")
    usable = Rect(mx, my, width - mx, height - my)
    nw, nh = normal
    yw, yh = ybranch
    pair_w, pair_h = nw + gx + yw, max(nh, yh)
    if pair_w > usable.width or pair_h > usable.height:
        raise ValueError("A corner pair does not fit inside the margins")
    corners = []
    for tag, right, top in [("BL", False, False), ("BR", True, False),
                            ("TL", False, True), ("TR", True, True)]:
        x = usable.right - pair_w if right else usable.left
        y = usable.top - pair_h if top else usable.bottom
        for kind, xx, w, h in [("normal", x, nw, nh), ("ybranch", x + nw + gx, yw, yh)]:
            yy = y + pair_h - h if top else y
            r = Rect(xx, yy, xx + w, yy + h)
            if any(not separated(r, other[2], gx, gy) for other in corners):
                raise ValueError("Corner references overlap or violate the minimum gap")
            corners.append((tag, kind, r))
    fw, fh = footprint
    nx = max(0, (usable.width + gx) // (fw + gx))
    ny = max(0, (usable.height + gy) // (fh + gy))
    available, excluded = [], []
    for row in range(ny):
        for col in range(nx):
            x, y = usable.left + col * (fw + gx), usable.bottom + row * (fh + gy)
            r = Rect(x, y, x + fw, y + fh)
            (available if all(separated(r, c[2], gx, gy) for c in corners)
             else excluded).append(r)
    return dict(usable=usable, nx=nx, ny=ny, corners=corners,
                available=available, excluded=excluded, footprint=footprint,
                pitch=(fw + gx, fh + gy))


def setup_layout(name):
    import pya
    import SiEPIC  # installs find_pin and other layout extensions
    from SiEPIC.scripts import load_klayout_library
    from SiEPIC.utils.layout import new_layout
    if not (PDK_ROOT / (tech_name + ".lyt")).exists():
        raise FileNotFoundError(f"Set PDK_ROOT to the installed ZEP technology folder: {PDK_ROOT}")
    if not pya.Technology.has_technology(tech_name):
        tech = pya.Technology.create_technology(tech_name)
        tech.load(str(PDK_ROOT / (tech_name + ".lyt")))
    # Load the updated Universal_GC together with native ZEP routing cells.
    for library, folder, pcells in [
        (tech_name, "SiEPICfab_EBeam_ZEP_fixed", "SiEPICfab_EBeam_ZEP_pcells"),
        (tech_name + "_Beta", "SiEPICfab_EBeam_ZEP_beta_fixed", "SiEPICfab_EBeam_ZEP_beta_pcells"),
    ]:
        if pya.Library.library_by_name(library, tech_name) is None:
            load_klayout_library(tech_name, library,
                                folder_gds="pymacros/" + folder,
                                folder_pcell="pymacros/" + pcells if pcells else None,
                                verbose=False)
    return new_layout(tech_name, name, GUI=False, overwrite=False)


def build_device(ly, name, settings=None, ybranch=False, label=None):
    """Reuse native cells, SiEPIC routing and the reference GC anchor convention."""
    import pya
    from SiEPIC.utils import create_cell2
    from SiEPIC.scripts import connect_pins_with_waveguide
    settings = dict(GC_SETTINGS if settings is None else settings)
    parameters = settings.get("parameters", {})
    if parameters:
        gc = ly.create_cell(settings["gc_cell_name"], settings["gc_library"], parameters)
        if gc is None or gc.bbox().empty():
            raise RuntimeError("Universal_GC is missing or empty; reload the updated ZEP PDK")
    else:
        gc = create_cell2(ly, settings["gc_cell_name"], settings["gc_library"])
    # The installed standard GC has a duplicate, identical PinRec path. Repair
    # only a local copy's metadata; fabrication shapes and library stay intact.
    clean_name = f"GC_routing_copy_{gc.cell_index()}"
    clean = ly.cell(clean_name)
    if clean is None:
        clean = ly.create_cell(clean_name)
        clean.copy_tree(gc)
        clean.flatten(True)
        seen = set()
        for shape in list(clean.shapes(ly.layer(ly.TECHNOLOGY["PinRec"])).each()):
            key = str(shape)
            if key in seen:
                shape.delete()
            else:
                seen.add(key)
    gc = clean
    gc_port = gc.find_pins()[0][0].pin_name
    cell = ly.create_cell(name)
    dbu = ly.dbu
    pitch = settings["gc_pitch"]
    if pitch < 2 * wg_bend_radius:
        raise ValueError("gc_pitch must accommodate two loopback bends")
    def inst(component, x, y):
        return cell.insert(pya.CellInstArray(component.cell_index(),
                           pya.Trans(round(x / dbu), round(y / dbu))))
    def route(a, pa, b, pb, turtle=None, turtle_b=None):
        result = connect_pins_with_waveguide(
            a, pa, b, pb, turtle_A=turtle, turtle_B=turtle_b,
            waveguide_type=WG_SETTINGS["wg_type"], r=wg_bend_radius,
            error_min_bend_radius=True, relaxed_pinnames=False)
        if not result:
            raise RuntimeError(f"Routing failed in {name}: {pa} -> {pb}")
    pin = gc.find_pin(gc_port).center
    x = gc.bbox().width() * dbu
    first = inst(gc, x, -pin.y * dbu)
    if not ybranch:
        second = inst(gc, x, pitch - pin.y * dbu)
        route(first, gc_port, second, gc_port, [75, 90, pitch, -90, 75])
    else:
        branch = create_cell2(ly, "ebeam_ybranch_te1550", tech_name + "_Beta")
        bp = branch.find_pin("opt1").center
        # Three equally spaced GCs: the center is already at y=0.
        upper = inst(gc, x, pitch - pin.y * dbu)
        lower = inst(gc, x, -pitch - pin.y * dbu)
        fp = first.find_pin(gc_port).center
        branch_inst = inst(branch, (fp.x - bp.x) * dbu + 30,
                           (fp.y - bp.y) * dbu)
        route(first, gc_port, branch_inst, "opt1")
        upper_y = upper.find_pin(gc_port).center.y * dbu
        lower_y = lower.find_pin(gc_port).center.y * dbu
        route(branch_inst, "opt2", upper, gc_port,
              [75, 90, upper_y - branch_inst.find_pin("opt2").center.y * dbu])
        route(branch_inst, "opt3", lower, gc_port,
              [75, -90, branch_inst.find_pin("opt3").center.y * dbu - lower_y])
        # GC envelopes must not intersect even within an alignment device.
        gcs = [first, upper, lower]
        for i, gc_instance in enumerate(gcs):
            for other in gcs[:i]:
                if not separated(rect_of(gc_instance), rect_of(other), 0, 0):
                    raise RuntimeError(f"GCs overlap inside alignment device {name}")
    text = pya.Text(label or f"opt_in_TE_1550_device_JS_{name}",
                    pya.Trans(round(x / dbu), 0))
    cell.shapes(ly.layer(ly.TECHNOLOGY["Text"])).insert(text).text_size = 10 / dbu
    if cell.bbox().empty():
        raise RuntimeError(f"Empty component: {name}")
    return cell


def size(cell):
    b = cell.bbox()
    return b.width(), b.height()


def make_plan(ly, normal, branch, footprint=None, config=None):
    cfg = dict(area_width=area_width, area_height=area_height,
               margin_x=margin_x, margin_y=margin_y,
               device_gap_x=device_gap_x, device_gap_y=device_gap_y,
               required_x_size=required_x_size, required_y_size=required_y_size)
    cfg.update(config or {})
    def units(v):
        if not math.isfinite(v) or v < 0:
            raise ValueError("Dimensions must be finite and nonnegative")
        return math.ceil(v / ly.dbu - 1e-9)
    fw, fh = footprint or size(normal)
    fw = max(fw, units(cfg["required_x_size"] or 0))
    fh = max(fh, units(cfg["required_y_size"] or 0))
    plan = packing(*(units(cfg[k]) for k in ["area_width", "area_height", "margin_x",
                   "margin_y", "device_gap_x", "device_gap_y"]),
                   (fw, fh), size(normal), size(branch))
    plan["config"] = cfg
    return plan


def place(top, component, rect):
    import pya
    b = component.bbox()
    if b.width() > rect.width or b.height() > rect.height:
        raise ValueError(f"{component.name} exceeds its allocated position")
    result = top.insert(pya.CellInstArray(component.cell_index(),
                        pya.Trans(rect.left - b.left, rect.bottom - b.bottom)))
    return result


def draw_plan(top, ly, plan, normal, branch, preview=True):
    import pya
    cfg = plan["config"]
    # Preview overlays are on explicitly named non-fabrication layers.
    boundary = ly.layer(ly.TECHNOLOGY["FloorPlan"])
    top.shapes(boundary).insert(pya.Box(0, 0, round(cfg["area_width"] / ly.dbu),
                                       round(cfg["area_height"] / ly.dbu)))
    for tag, kind, r in plan["corners"]:
        # Distinct corner labels, with identical label sizing/anchor convention.
        source = normal if kind == "normal" else branch
        copy = ly.create_cell(f"ALIGN_{tag}_{kind}")
        copy.copy_tree(source)
        text_layer = ly.layer(ly.TECHNOLOGY["Text"])
        for shape in copy.shapes(text_layer).each():
            if shape.is_text():
                text = shape.text
                text.string = f"opt_in_TE_1550_device_JS_ALIGNMENT_{tag}_{kind}"
                shape.text = text
        # Text is metadata, but still included in the measured bounds.
        place(top, copy, r)
    if preview:
        for key, layer, label in [("available", 990, "SWEEP_AVAILABLE"),
                                   ("excluded", 991, "CORNER_EXCLUDED")]:
            li = ly.layer(pya.LayerInfo(layer, 0, label))
            for r in plan[key]:
                top.shapes(li).insert(pya.Box(r.left, r.bottom, r.right, r.top))


def report(ly, plan, normal):
    d = ly.dbu
    cfg = plan["config"]
    print(f"Fixed layout area: {cfg['area_width']} um x {cfg['area_height']} um")
    print(f"Usable area: {plan['usable'].width*d:g} um x {plan['usable'].height*d:g} um")
    print(f"GC loopback bounding box: {normal.bbox().to_dtype(d)} um")
    print(f"GC loopback width / height: {size(normal)[0]*d:g} / {size(normal)[1]*d:g} um")
    print(f"Allocated footprint: {plan['footprint'][0]*d:g} x {plan['footprint'][1]*d:g} um")
    print(f"Device pitch in X / Y: {plan['pitch'][0]*d:g} / {plan['pitch'][1]*d:g} um")
    print(f"Number fitting in X / Y: {plan['nx']} / {plan['ny']}")
    print(f"Standard grid capacity before corner reservation: {plan['nx']*plan['ny']}")
    print(f"Corner/reference devices: {len(plan['corners'])} (4 normal + 4 Y-branch)")
    print(f"Grid positions excluded by corner references: {len(plan['excluded'])}")
    print(f"Remaining sweep positions: {len(plan['available'])}")
    print(f"Total device capacity including references: {len(plan['available'])+8}")


def export(top, name, show=True, output_dir=None):
    from SiEPIC.scripts import export_layout
    from SiEPIC._globals import Python_Env
    from SiEPIC.utils import klive
    folder = Path(output_dir) if output_dir else BASE_DIR / "gds_in"
    folder.mkdir(parents=True, exist_ok=True)
    result = export_layout(top, str(folder), name, relative_path="",
                           format="gds", screenshot=False)
    if show and Python_Env == "Script":
        klive.show(result, technology=tech_name)
    return result


def rect_of(instance):
    b = instance.bbox()
    return Rect(b.left, b.bottom, b.right, b.top)


def validate_placements(top, ly, plan):
    """Check all real devices, including alignment-to-formal and corner pairs."""
    cfg = plan["config"]
    gx = math.ceil(cfg["device_gap_x"] / ly.dbu - 1e-9)
    gy = math.ceil(cfg["device_gap_y"] / ly.dbu - 1e-9)
    bounds = plan["usable"]
    placed = [(i.cell.name, rect_of(i)) for i in top.each_inst()]
    for index, (name, r) in enumerate(placed):
        if not (bounds.left <= r.left and r.right <= bounds.right
                and bounds.bottom <= r.bottom and r.top <= bounds.top):
            raise RuntimeError(f"{name} exceeds the usable area")
        for other_name, other in placed[:index]:
            if not separated(r, other, gx, gy):
                raise RuntimeError(f"Overlap/gap violation: {name} and {other_name}")
    print(f"PASS: {len(placed)} real devices inside margins; no overlaps or gap violations")
    return len(placed)


# %% Generate all legal normal loopbacks with consecutive one-based numbers.
def main(show=show_layout, output_dir=None):
    import csv
    top, ly = setup_layout("Si_GC_capacity")
    normal = build_device(ly, "GC_standard")
    branch = build_device(ly, "GC_Y_reference", ybranch=True)
    plan = make_plan(ly, normal, branch)
    report(ly, plan, normal)
    draw_plan(top, ly, plan, normal, branch, preview=False)
    rows = []
    # Share routed geometry; each wrapper has its own real cell name and label.
    text_layer = ly.layer(ly.TECHNOLOGY["Text"])
    normal.shapes(text_layer).clear()
    import pya
    for number, rect in enumerate(plan["available"], 1):
        name = f"GC_LOOPBACK_{number:04d}"
        cell = ly.create_cell(name)
        cell.insert(pya.CellInstArray(normal.cell_index(), pya.Trans()))
        label = pya.Text(f"opt_in_TE_1550_device_JS_{name}",
                         pya.Trans(round(normal.bbox().center().x), 0))
        cell.shapes(text_layer).insert(label).text_size = 10 / ly.dbu
        instance = place(top, cell, rect)
        b = instance.bbox()
        rows.append(dict(number=number, cell_name=name, x_um=b.left*ly.dbu,
                         y_um=b.bottom*ly.dbu, width_um=b.width()*ly.dbu,
                         height_um=b.height()*ly.dbu))
    count = validate_placements(top, ly, plan)
    assert count == len(rows) + 8
    folder = Path(output_dir) if output_dir else BASE_DIR / "gds_in"
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "Si_GC_capacity.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["number", "cell_name", "x_um",
                                "y_um", "width_um", "height_um"])
        writer.writeheader()
        writer.writerows(rows)
    plan["loopbacks"] = rows
    print(f"Drawn normal loopbacks: {len(rows)}; alignment devices: 8; total: {count}")
    if rows:
        print(f"Numbered cells: {rows[0]['cell_name']} through {rows[-1]['cell_name']}")
    export(top, "Si_GC_capacity", show, output_dir)
    return top, ly, plan


if __name__ == "__main__":
    top_cell, ly, capacity = main()
