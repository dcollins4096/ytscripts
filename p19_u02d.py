
import pyximport; pyximport.install()
import particle_ops
import clump_particles
reload(clump_particles)

if 1:
  thisdir = '/scratch1/dcollins/Paper19/B02/u02-128-d-alpha0.5'
  frame = 360
  width = (0.05,'code_length')
  setname = '%s/DD%04d/data%04d'%(thisdir,frame,frame)
  ds = yt.load(setname)
  width = ds.arr(0.05,'code_length')
  val,loc = ds.find_max('density')
  #loc = ds.arr([0.5]*3,'code_length')

  sphere = ds.sphere(loc,width)
  proj2 = ds.proj('density',0,data_source=sphere,center=loc)
  pw = proj2.to_pw(center = loc, width = (0.1,'code_length'))
  pw.save('test2.png')

  if 0:
    proj_full = ds.proj('density',0, center = loc) #width = (1.0,'code_length'))
    pw_full = proj_full.to_pw(center = loc,width=(1.0,'code_length'))
    pw_full.annotate_particles(1.0,stride=16)
    pw_full.save('test_full.png')

if 1:
  indices_late,xpos_late,ypos_late,zpos_late = clump_particles.particles_from_clump(sphere)

if 1:
  pw.annotate_dave_particles(1.0, indices=indices_late)
  pw.save('test3.png')





