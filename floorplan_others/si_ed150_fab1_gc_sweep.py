# %%
"""Universal GC parameter sweep using si_ed150_fab1_gc_capacity.py.

The shared module supplies routing, measured packing, eight corner alignment
structures (including three-GC Y branches), gap validation and GDS export.
All requested variants are measured before choosing a common grid footprint.
Unused positions remain empty; excess combinations are reported and omitted.
Run as Python or VS Code Jupyter cells. Geometry dimensions are micrometers.
"""
from pathlib import Path
import importlib.util
from itertools import product
from copy import deepcopy
import csv
import sys

BASE_DIR = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
# Support Interactive Window execution from either the repo or notebook folder.
module_path = BASE_DIR / "si_ed150_fab1_gc_capacity.py"
if not module_path.exists():
    module_path = BASE_DIR / "floorplan_BTOED150_tapeout0" / module_path.name
spec = importlib.util.spec_from_file_location("si_gc_capacity", module_path)
capacity = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = capacity
spec.loader.exec_module(capacity)
BASE_DIR = module_path.parent

# %% Layout settings inherit the capacity defaults, including the reference GCs.
area_width = capacity.area_width
area_height = capacity.area_height
margin_x = capacity.margin_x
margin_y = capacity.margin_y
device_gap_x = capacity.device_gap_x
device_gap_y = capacity.device_gap_y
required_x_size = capacity.required_x_size
required_y_size = capacity.required_y_size
GC_SETTINGS = deepcopy(capacity.GC_SETTINGS)
# Corner/alignment GCs: supplied pitches are nm, stored here in um.
ALIGNMENT_GC_SETTINGS = deepcopy(GC_SETTINGS)
ALIGNMENT_GC_SETTINGS["parameters"].update(
    target_length=20.0, duty_cycle=0.70, pitch=0.755,
    radius=15.0, y_span=10.0, L_extra=5.0,
)
main_sweep_parameters = dict(
    target_length=[15.0, 20.0, 30.0],
    duty_cycle=[value / 100 for value in range(40, 76, 5)],
    pitch=[value / 1000 for value in range(675, 906, 10)],
    radius=[15.0, 25.0], y_span=[10.0, 15.0], L_extra=[5.0],
)
# Each additional group varies only one parameter from the alignment defaults.
additional_sweeps = [
    ("L_extra", [0.0, 3.0, 5.0, 7.0], 2),
    ("y_span", [7.5, 10.0, 12.5, 15.0], 2),
    ("radius", [10.0, 15.0, 20.0, 25.0, 30.0, 35.0], 2),
    ("target_length", [15.0, 20.0, 25.0, 30.0], 2),
]
show_layout = True


def combinations():
    keys = list(main_sweep_parameters)
    for values in product(*(main_sweep_parameters[key] for key in keys)):
        settings = deepcopy(ALIGNMENT_GC_SETTINGS)
        settings["parameters"].update(dict(zip(keys, values)))
        settings.update(sweep_group="main", repeat=1)
        yield settings
    for parameter, values, repeats in additional_sweeps:
        for repeat in range(1, repeats + 1):
            for value in values:
                settings = deepcopy(ALIGNMENT_GC_SETTINGS)
                settings["parameters"][parameter] = value
                settings.update(sweep_group=parameter, repeat=repeat)
                yield settings


# %% Measure real routed variants, reserve alignment cells, then fill legal slots.
def main(show=show_layout, output_dir=None):
    top, ly = capacity.setup_layout("Si_GC_loopback_sweep")
    normal = capacity.build_device(ly, "GC_standard", ALIGNMENT_GC_SETTINGS)
    branch = capacity.build_device(ly, "GC_Y_reference", ALIGNMENT_GC_SETTINGS, ybranch=True)
    requested = []
    settings_list = list(combinations())
    group_counts = {}
    for settings in settings_list:
        group = settings["sweep_group"]
        group_counts[group] = group_counts.get(group, 0) + 1
    print(f"Requested groups: {group_counts}; total: {len(settings_list)}")
    for number, settings in enumerate(settings_list, 1):
        params = settings["parameters"]
        name = (
            f"GC{number}_Len{params['target_length']:g}"
            f"DC{params['duty_cycle']*100:g}pitch{params['pitch']:.3f}"
            f"Rad{params['radius']:g}Ys{params['y_span']:g}Lex{params['L_extra']:g}"
        ).replace(".", "p")
        label = f"opt_in_TE_1550_device_JS_GCsweep_{name}"
        device = capacity.build_device(ly, name, settings, label=label)
        requested.append((number, device, settings))
    measured = [normal] + [item[1] for item in requested]
    footprint = tuple(max(capacity.size(cell)[axis] for cell in measured)
                      for axis in (0, 1))
    cfg = {key: globals()[key] for key in ["area_width", "area_height", "margin_x",
           "margin_y", "device_gap_x", "device_gap_y", "required_x_size", "required_y_size"]}
    plan = capacity.make_plan(ly, normal, branch, footprint, cfg)
    capacity.report(ly, plan, normal)
    capacity.draw_plan(top, ly, plan, normal, branch, preview=False)
    available = len(plan["available"])
    print(f"Requested sweep structures: {len(requested)}; available positions: {available}")
    if len(requested) > available:
        print(f"WARNING: {len(requested)-available} requested structures cannot be placed")
    rows = []
    for (number, device, settings), rect in zip(requested, plan["available"]):
        instance = capacity.place(top, device, rect)
        b = instance.bbox()
        rows.append(dict(number=number, device_cell=device.name,
                         sweep_group=settings["sweep_group"], repeat=settings["repeat"],
                        #  gc_cell_name=settings["gc_cell_name"],
                        #  gc_library=settings["gc_library"], 
                         gc_pitch=settings["gc_pitch"],
                         **settings["parameters"], x_um=b.left*ly.dbu,
                         y_um=b.bottom*ly.dbu, width_um=b.width()*ly.dbu,
                         height_um=b.height()*ly.dbu))
    count = capacity.validate_placements(top, ly, plan)
    if count != len(rows) + 8:
        raise RuntimeError("Unexpected number of placed sweep/alignment devices")
    folder = Path(output_dir) if output_dir else BASE_DIR / "designs"
    folder.mkdir(parents=True, exist_ok=True)
    fields = ["number", "device_cell", "sweep_group", "repeat", "gc_cell_name", "gc_library", "gc_pitch",
              *GC_SETTINGS["parameters"], "x_um", "y_um", "width_um", "height_um"]
    with (folder / "Si_GC_loopback_sweep.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    placed_numbers = {row["number"] for row in rows}
    with (folder / "Si_GC_loopback_sweep_requested.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["number", "device_cell", "sweep_group",
            "repeat", "placed", *GC_SETTINGS["parameters"]])
        writer.writeheader()
        for number, device, settings in requested:
            writer.writerow(dict(number=number, device_cell=device.name,
                sweep_group=settings["sweep_group"], repeat=settings["repeat"],
                placed=number in placed_numbers, **settings["parameters"]))
    for group in group_counts:
        placed_count = sum(row["sweep_group"] == group for row in rows)
        print(f"{group}: requested {group_counts[group]}, placed {placed_count}, "
              f"unplaced {group_counts[group]-placed_count}")
    capacity.export(top, "Si_GC_loopback_sweep", show, folder)
    print(f"Placed: {len(rows)} sweep + 8 alignment devices; empty positions: {available-len(rows)}")
    return top, ly, plan, rows


if __name__ == "__main__":
    top_cell, ly, packing_plan, sweep_rows = main()
