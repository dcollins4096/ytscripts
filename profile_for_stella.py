
        prof_tracer = yt.create_profile(ds.all_data(),('deposit','all_cic'),fields='cell_volume',weight_field=None)
        plt.plot(0.5*(prof_tracer.x_bins[1:]+prof_tracer.x_bins[0:-1])*scaler,prof_tracer['cell_volume'],c='r',
                 label = 'tracer')
