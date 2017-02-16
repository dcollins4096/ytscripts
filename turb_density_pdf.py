execfile('go')
car = taxi.taxi(sys.argv[1])
car.profile(['density','cell_volume'])
