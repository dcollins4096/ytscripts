if 'ef' not in dir():
    execfile('go')
reload(taxi)
#for name in ['az20','az21','az22']: #['aa19']: #,'aa20','aa21','aa22']:
#    car = taxi.taxi(name)
#    car.qb_load()
#    car.qb.plot()
if 'all_proj' not in dir():
    all_proj = open('all_proj','r')

    all_files = []
    for line in all_proj:
        all_files.append(line.split("/")[-1][:-1])

counter = 0

for car in flt:
    fptr = open("plot_log_%s.txt"%car.outname,"a+")
    car.frames='every 10'
    for frame in car.return_frames():
        for field in car.fields:
            for axis in [0,1,0]:
                counter += 1
                if counter <0:
                    continue
                #ax22_0010_Projection_z_magnetic_field_strength.png
                fname = "%s_%04d_Projection_%s_%s.png"%(car.name,frame,'xyz'[axis],field)
                print fname in all_files, fname
                if fname in all_files:
                    fptr.write("%s %s %s %s\n"%(str(frame),field,str(axis),'Full'))
    fptr.close()
            
if 0:
    car.outname = 'clobber_test'
    car.operation='CenterSlice'
    car.frames=[600]
    car.axis=0
    car.fields=['density','Bx']

    car.clobber=True
    car.plot()
