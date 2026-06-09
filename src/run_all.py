"""
run_all.py : regenerate every figure and GIF in the repository.

    cd src && python3 run_all.py

Runs the three models in order and writes all outputs to ../figures/.
"""
import importlib

for mod in ("reaction_diffusion", "hedgehog_assembly", "corncob_coaggregation",
            "field_3d", "reaction_diffusion_3d", "succession_3d"):
    print(f"\n=== running {mod} ===")
    m = importlib.import_module(mod)
    # each module guards its work under __main__, so call the work explicitly
    if mod == "reaction_diffusion":
        snaps = m.simulate(); m.make_gif(snaps); m.make_static(snaps)
    elif mod == "hedgehog_assembly":
        h, _ = m.build(); m.radial_profile(h)
    elif mod == "corncob_coaggregation":
        cc, _ = m.build(); m.analysis(cc)
    elif mod == "field_3d":
        m.single_gif("O$_2$", "D_field_O2_3d.gif")
        m.single_gif("CO$_2$", "D_field_CO2_3d.gif")
        m.single_gif("lactate", "D_field_lactate_3d.gif")
        m.multipanel_gif()
    elif mod == "reaction_diffusion_3d":
        cells, fils = m.build_biomass()
        Wresp = m.deposit(cells, m.RESP_W); Wferm = m.deposit(cells, m.FERM_W)
        O, Lf, Cf = m.solve(Wresp, Wferm)
        m.gif_biomass(cells, fils, "E_biomass_3d.gif")
        m.gif_field(O/max(O.max(), 1e-9), "turbo", "O$_2$ field (3D, coupled)", "E_field3d_O2.gif", gamma=1.9)
        m.gif_field(Cf, "viridis", "CO$_2$ field (3D, coupled)", "E_field3d_CO2.gif", gamma=1.5)
        m.gif_field(Lf, "YlGn", "lactate field (3D, coupled)", "E_field3d_lactate.gif", gamma=1.3)
        m.slices_png(O, cells); m.heterogeneity_png(O); m.multipanel(cells, fils, O)
        m.topview_png(cells, fils, O)
    else:  # succession_3d
        (pos, tax, tapp), fils, snaps, metrics = m.run_succession()
        m.gif_topview(pos, tax, tapp, fils, snaps)
        m.gif_multipanel(pos, tax, tapp, fils, snaps, metrics)
        m.gif_species(pos, tax, tapp, fils, snaps)
        m.gif_single_field(snaps, "O", m.CMAP["O2"], "O$_2$ (fraction of saliva)", "F_succession_O2_3d.gif")
        m.gif_single_field(snaps, "C", m.CMAP["CO2"], "CO$_2$ (normalised)", "F_succession_CO2_3d.gif")
        m.gif_single_field(snaps, "L", m.CMAP["lactate"], "lactate (normalised)", "F_succession_lactate_3d.gif")
        m.plot_dynamics(metrics)

print("\n=== running analysis dashboard ===")
import analysis
analysis.main()

print("\nAll figures written to ../figures/")
