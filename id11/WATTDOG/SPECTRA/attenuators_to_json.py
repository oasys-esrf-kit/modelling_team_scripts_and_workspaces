
import json

density = { 'C':  3.52,
            'Cu': 8.96,
            'Au': 19.32,
            'Ag': 10.49,
            'Al': 2.70,
            'Si': 2.33 }

EMPTY = {
    "name": "Empty",
    "substance": "Air, Dry (near sea level)",
    "thickness": 0.0,
    "density": 1.0e-6,
}

def foil_entry(substance, thickness_mm):
    return {
        "name": "{} filter {:.4g} mm".format(substance, thickness_mm),
        "substance": substance,
        "thickness": thickness_mm,
        "density": density[substance],
    }

def axes_to_json(axes, add_empty=False):
    doc = {}
    for i, axis in enumerate(axes):
        axis_key = "Axis{}".format(i + 1)
        axis_doc = {"_name": "Att{}".format(i + 1), "_att_pos": 0}
        offset = 1
        if add_empty:
            axis_doc["filter1"] = EMPTY
            offset = 2
        for j, (substance, thickness) in enumerate(axis):
            axis_doc["filter{}".format(j + offset)] = foil_entry(substance, thickness)
        doc[axis_key] = axis_doc
    return doc

# extended syned
from orangecontrib.esrf.syned.util.syned_filter_with_density import FilterWithDensity
from orangecontrib.esrf.syned.util.syned_filter_packs import FilterBox, FilterBlock
def _att_dic_to_syned_filterbox(att_dic):
    n_keys = 0
    keys = []
    for key in att_dic.keys():
        n_keys += 1
        keys.append(key)

    # update combo boxes
    block_list = []
    for i in range(n_keys):
        items = []
        for filter in att_dic[keys[i]].keys():
            if filter[0] != "_":
                item = att_dic[keys[i]][filter]
                f = FilterWithDensity(name=item['name'],
                                      material=item['substance'],
                                      thickness=item['thickness'],
                                      density=item['density'])
                items.append(f)

        block_list.append(FilterBlock(filters_list=items))

    return FilterBox(filter_blocks_list=block_list)

# ---------- 2026 ----------
cvd = ( 'C', 0.3 )
D   = cvd,
Cu  = ('Cu', 0.0025), cvd, ('Cu', 0.0025)
d   = ('Au', 0.0038), cvd, ('Au', 0.0038)
c   = ('Au', 0.0013), cvd, ('Au', 0.0038)
Al1 = ('Al', 1.0),

axes_2026 = [D, Cu, Cu, d, c, c, Al1, Al1]

print("***** 2026 *****")
for i, ax in enumerate(axes_2026):
    print("axis ", i+1, ": ", ax)

out_2026 = "id11_wattdog_attenuators_2026.json"
with open(out_2026, 'w') as f:
    json.dump(axes_to_json(axes_2026), f, indent=4)
print("Written", out_2026)


# ---------- 2028 (Post Refurbishment) ----------
cvd5 = ('C', 0.5)

ax1 = cvd5, ('Cu', 0.0015), cvd, ('Cu', 0.0015)
ax2 = ('Cu', 0.0025), cvd, ('Cu', 0.0025), ('Cu', 0.005),  cvd, ('Cu', 0.005)
ax3 = ('Cu', 0.0075), cvd, ('Cu', 0.0075), ('Cu', 0.010),  cvd, ('Cu', 0.010)
ax4 = ('Cu', 0.0150), cvd, ('Cu', 0.0150), ('Ag', 0.005),  cvd, ('Ag', 0.005)
ax5 = ('Cu', 0.0075), cvd, ('Cu', 0.0075), ('Ag', 0.010),  cvd, ('Ag', 0.010)
ax6 = ('Cu', 0.0150), cvd, ('Cu', 0.0150), ('Ag', 0.025),  cvd, ('Ag', 0.025)

axes_2028 = [ax1, ax2, ax3, ax4, ax5, ax6]

print("***** 2028 *****")
for i, ax in enumerate(axes_2028):
    print("axis ", i+1, ": ", ax)

out_2028 = "id11_wattdog_attenuators_2028.json"
with open(out_2028, 'w') as f:
    json.dump(axes_to_json(axes_2028), f, indent=4)
print("Written", out_2028)


#
# read files to dictionary
#
is_remote = 0
syned_file_name = out_2028
import urllib
if is_remote:
    response = urllib.request.urlopen(syned_file_name)
    att_dic = json.load(response)
else:
    with open(syned_file_name) as att_file:
        att_dic = json.load(att_file)

print(att_dic)

#
# from json dictionary to syned opject
#
syned_filterbox = _att_dic_to_syned_filterbox(att_dic)
# print(syned_filterbox.info())
syned_filterbox.to_json("tmp.json")


#
# read syned json files
#
from syned.util.json_tools import load_from_json_file, load_from_json_url
is_remote = 0
syned_file_name = "tmp.json"
_exec = ("from orangecontrib.esrf.syned.util.syned_filter_with_density import FilterWithDensity\n"
         "from orangecontrib.esrf.syned.util.syned_filter_packs import FilterBlock, FilterBox")
if is_remote:
    syned_filterbox2 = load_from_json_url(syned_file_name, exec_commands=_exec)
else:
    syned_filterbox2 = load_from_json_file(syned_file_name, exec_commands=_exec)

m, t, d = syned_filterbox2.get_lists_materials_thicknesses_densities(cumulate=1)
for i in range(len(m)):
    print(i, m[i], t[i], d[i])
# print("materials: ", m)
# print("thicknesses: ", t)
# print("densities: ", d)