
if 'ef' not in dir():
    execfile('go')
reload(taxi)
if 'always' not in dir():
    always=True
if 'garr' not in dir() or always:
    r15=taxi.taxi('r15')
    r15x=taxi.taxi('r15b')
    r17b=taxi.taxi('r17b')
    r17=taxi.taxi('r17')
    flt = taxi.fleet([r15x,r17b])#,r15,r17])
    flt('car.load(0)')
    garr = flt('output.append( car.ds.index.grids[0])')
for a,b in [['d','Density'],['vx','x-velocity'],['vy','y-velocity'],['vz','z-velocity'],['E','TotalEnergy']]:
    quan[b]=quan[a]
def moo(field):
    print( "%5s            %0.16f \pm %0.16e "%(flt[0].name, np.mean(garr[0][field].v), np.std(garr[0][field].v)))
    print( "%5s            %0.16f \pm %0.16e "%(flt[1].name, np.mean(garr[1][field].v), np.std(garr[1][field].v)))
    #print( "%5s            %0.16f \pm %0.16e "%(flt[3].name, np.mean(garr[3][field].v), np.std(garr[3][field].v)))
    print( "%5s %5s diff %0.16e"%(flt[0].name, flt[1].name,(np.abs(garr[0][field].v-garr[1][field].v)).sum()))
    #print( "%5s %5s diff %0.16e"%(flt[3].name, flt[1].name,(np.abs(garr[3][field].v-garr[1][field].v)).sum()))
    #print( "%5s %5s diff %0.16e"%(flt[3].name, flt[0].name,(np.abs(garr[3][field].v-garr[0][field].v)).sum()))
    #print( "2 q g2 %0.16f \pm %0.16e quan %0.16e diff %0.16e"%( np.mean(garr[2][field].v), np.std(np.abs(garr[2][field].v)), quan[field], (np.abs(garr[2][field].v-quan[field])).sum()))

for field in ['Density','x-velocity','y-velocity','z-velocity','TotalEnergy']:
    print("==== %s ===="%field)
    moo(field)
