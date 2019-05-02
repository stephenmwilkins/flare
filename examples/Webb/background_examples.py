from jwst_backgrounds import jbt
 
pointings = [{'ra':53.1625,'dec':27.7914,'name':'HUDF'},
             {'ra':83.8187,'dec':-5.3872,'name':'Trapezium'},
             {'ra':266.41683,'dec':-29.00781,'name':'Galactic Center'}]
 
wavelength = 4.4 #micron
 
for pointing in pointings:
    bg = jbt.background(pointing['ra'],pointing['dec'],wavelength)
    print(pointing['name'], bg.bathtub['total_thiswave'][0], 'MJy/sr')
    
    
    
jbt.get_background(223.555, -54.395, 2.15, thresh=1.1, plot_background=True, plot_bathtub=True, write_bathtub=True)