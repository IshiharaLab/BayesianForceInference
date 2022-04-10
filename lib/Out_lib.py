import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import   Polygon
from matplotlib.collections import PatchCollection

'''
   Functions for output rsults and for drawing cells and reults
     DrawCells(x,y,edge,cell):
     Draw_Tension(x,y,T,T_LINE_WIDTH = 2.0, savefile = [])
     Draw_Pressure(x,y,P,savefile = []):
'''


#######
def OutputresultsTP(filename,outfile,x,y,edge,cell,T,P,mu):
    oFid = open( outfile, 'w' )
    oFid.write('### results saved in %s \n' % (filename) )
    #ofid.wrute( '### ABIC= %e ' % ( ABIC ) )
    oFid.write("### mu=   %e \n" % ( mu ) )

    oFid.write("# edge tension : edge_id, inf.tension, - -, (x1 y2), (x2 y2) \n")
    for (i, e) in enumerate(edge):
        oFid.write("%3d   %e - -    (%.3f %.3f) (%.3f %.3f) \n" % (i,T[i],e.x1,e.y1, e.x2,e.y2) )

    oFid.write("\n")
    oFid.write("# cell pressure : - - cell_id, inf.pressure\n")
    for (i, c) in enumerate(cell):
        oFid.write("- - %3d   %e\n" % (i,P[i]) )

    oFid.write("\n")
    oFid.close()

######
def DrawCells(x,y,edge,cell):
    print('...    Show cells               ')

    fP, ax = plt.subplots(1,figsize=(10,10))
    fP.patch.set_facecolor('white')

    for (i,c) in enumerate(cell):
        pg = np.vstack( (x[c.junc],y[c.junc]) ).T
        #pg = np.vstack( (cinf[i][0::2],cinf[i][1::2]) ).T
        if(c.in_out=='i'):
            polygon = Polygon( pg, closed = True, fc = 'gold' )
            ax.add_patch(polygon)

    for e in edge:
        plt.plot([e.x1,e.x2],[e.y1,e.y2],'-k',linewidth=1)

    ax.autoscale_view()
    ax.set_aspect('equal', 'datalim')
    plt.axis('off')

    #plt.pause(0.2)
    #plt.savefig("tmpCell.eps")



##################
def Draw_Tension(x,y,T,edge,T_LINE_WIDTH = 2.0, tmin = 0,tmax=2.0, savefile = ''):
    print('  ...   Show tensions, saved as ',savefile)
    #plt.rcParams['figure.figsize'] = 16,16
    fT, ax = plt.subplots(figsize=(16,16))
    fT.patch.set_facecolor('white')
    # maxT = max(T)
    # minT = min(T)
    lines = []
    colors = []

    for (i,e) in enumerate(edge):
        lines.append( [(e.x1,e.y1), (e.x2,e.y2)] )
        colors.append( T[i] )

    line_segments = LineCollection( lines, linewidths=T_LINE_WIDTH, linestyles='solid', cmap = matplotlib.cm.jet)
    line_segments.set_array( np.array(colors) )
    line_segments.set_clim([tmin, tmax])
    ax.add_collection(line_segments)
    ax.autoscale_view()
    ax.set_aspect('equal', 'datalim')
    plt.axis('off')
    fT.colorbar(line_segments, ax=ax, orientation='horizontal')

    fT.show()
    if savefile != '':
        fT.savefig(savefile,format = 'eps')



def Draw_Pressure(x,y,P,edge,cell,pmin ,pmax,savefile = ''):
    print('  ...   Show pressures, saved as ',savefile)
    fP, ax = plt.subplots(1,figsize=(16,16))
    fP.patch.set_facecolor('white')
    #maxP = max(P);
    #minP = min(P);
    patches = []
    colors = []

    for (i,e) in enumerate(edge):
        plt.plot([e.x1,e.x2],[e.y1,e.y2],linewidth=1,color = "0.2")

    for (i,c) in enumerate(cell):
        if  c.in_out == 'i':
            pg = np.vstack( (x[c.junc], y[c.junc]) ).T
            polygon = Polygon( pg, closed = True )
            patches.append(polygon)
            colors.append(P[i])

    collection = PatchCollection(patches, linewidths=0.5, cmap = matplotlib.cm.cool)
    collection.set_array( np.array(colors) )
    collection.set_clim([pmin, pmax])
    ax.add_collection( collection )
    ax.autoscale_view()
    ax.set_aspect('equal', 'datalim')
    plt.axis('off')

    fP.colorbar(collection, ax=ax,orientation='horizontal')
    fP.show()
    if savefile != '':
        fP.savefig( savefile, format = 'eps' )

