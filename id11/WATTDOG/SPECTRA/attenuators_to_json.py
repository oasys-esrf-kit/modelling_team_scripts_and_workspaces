import json, os
from syned.beamline.optical_elements.absorbers.filter_box import FilterBox

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

if __name__ in ["__main__", "builtins"]:

    # wattdog ideas: https://confluence.esrf.fr/pages/viewpage.action?pageId=176196123

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



    """
    https://confluence.esrf.fr/display/PROJID11REFURB/TDR+documents?preview=%2F180945561%2F180945563%2FID11-Refurbishment-TDR_V5.pdf
    ID11 currently has  six attenuator axes that can each host 2 CVD diamond disks coated with metals on 
    either side. This allows up to 12 disks to be inserted in the beam, with a wide range of combinations possible.
    The scheme, which is typically in normal use at ID11 for the Laue mono  today, is:
    
    1 = D = CVD diamond, 300 um,
    2 = Cu = CVD diamond, 300 um, with 5 um copper (2.5 um on each side)
    3 = Cu = CVD diamond, 300 um, with 5 um copper (2.5 um on each side)
    4 = d = CVD diamond, 300 um, with 7.6 um gold (3.8 um on each side)
    5 = c = CVD diamond, 300 um, with 5.1 um gold (1.3 upstream, 3.8 um downstream)
    6 = c = CVD diamond, 300 um, with 5.1 um gold (1.3 upstream, 3.8 um downstream)
    7 = Aluminium 1 mm
    8 = Aluminium 2 mm
    9 = Aluminium 4 mm
    
    """

    # ---------- 2026 ---------- https://confluence.esrf.fr/display/ID11KB/04+-+Attenuator

    cvd = ( 'C', 0.3 )

    FE_window = cvd,
    D   = cvd,
    Cu  = ('Cu', 0.0025), cvd, ('Cu', 0.0025)
    d   = ('Au', 0.0038), cvd, ('Au', 0.0038)
    c   = ('Au', 0.0013), cvd, ('Au', 0.0038)
    Al1 = ('Al', 1.0),

    axes_2026 = [FE_window, D, Cu, Cu, d, c, c, Al1, Al1]

    print("***** 2026 *****")
    for i, ax in enumerate(axes_2026):
        print("axis ", i+1, ": ", ax)

    out_2026 = "id11_wattdog_attenuators_2026.json"
    with open(out_2026, 'w') as f:
        json.dump(axes_to_json(axes_2026), f, indent=4)
    print("Written", out_2026)


    # ---------- 2028 (Post Refurbishment) ----------   https://confluence.esrf.fr/display/PROJID11REFURB/Attenuators

    ax0 = ('C', 0.3), # window at the front end
    ax1 = ('C', 0.5), ('Cu', 0.0015), cvd, ('Cu', 0.0015)
    ax2 = ('Cu', 0.0025), cvd, ('Cu', 0.0025), ('Cu', 0.005),  cvd, ('Cu', 0.005)
    ax3 = ('Cu', 0.0075), cvd, ('Cu', 0.0075), ('Cu', 0.010),  cvd, ('Cu', 0.010)
    ax4 = ('Cu', 0.0150), cvd, ('Cu', 0.0150), ('Ag', 0.005),  cvd, ('Ag', 0.005)
    ax5 = ('Cu', 0.0075), cvd, ('Cu', 0.0075), ('Ag', 0.010),  cvd, ('Ag', 0.010)
    ax6 = ('Cu', 0.0150), cvd, ('Cu', 0.0150), ('Ag', 0.025),  cvd, ('Ag', 0.025)

    axes_2028 = [ax0, ax1, ax2, ax3, ax4, ax5, ax6]

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
    check_output = 1
    if check_output:
        file_name = out_2028
        is_remote = 0
        for file_name in [out_2026, out_2028]:
            if is_remote:
                import urllib
                response = urllib.request.urlopen(file_name)
                att_dic = json.load(response)
            else:
                with open(file_name) as att_file:
                    att_dic = json.load(att_file)

            # print(att_dic)

            #
            # from json dictionary to syned object
            #
            syned_filterbox = FilterBox.from_plane_json_dict(att_dic)
            #
            # write json syned
            #
            syned_file_name = os.path.splitext(file_name)[0]+"_syned.json"
            syned_filterbox.to_json(syned_file_name)

            #
            # read json syned
            #
            from syned.util.json_tools import load_from_json_file, load_from_json_url
            if is_remote:
                syned_filterbox2 = load_from_json_url(syned_file_name)
            else:
                syned_filterbox2 = load_from_json_file(syned_file_name)

            m, t, d = syned_filterbox2.get_lists_materials_thicknesses_densities(cumulate=1)
            print("************* cumulated:  ", syned_file_name)
            for i in range(len(m)):
                print(i, m[i], t[i], d[i])
