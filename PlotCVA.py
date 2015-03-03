#!/usr/bin/python
''' 
Plotting the CVA in 2D and 3D
'''
# importing bit##########################################################################################
import matplotlib, sys
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse
from numpy import mean,var,sqrt
from mpl_toolkits.mplot3d import Axes3D
import cPickle as pickle
from scipy import stats
from random import shuffle
# ########################################################################################################
colors=[
    '#ED9121', '#EE8262', '#EE1289', '#556B2F', '#FF8C00', '#8B7B8B', '#0000EE', '#EED5D2', 
    '#BA55D3', '#912CEE', '#2F4F4F', '#D15FEE', '#008B8B', '#B23AEE', '#8B7765', '#54FF9F',
    '#8B8386', '#FF4040', '#EEA9B8', '#388E8E', '#6E8B3D', '#33A1C9', '#EE3A8C', '#FF00FF',
    '#436EEE', '#8B864E', '#808000', '#1874CD', '#BCD2EE', '#A9A9A9', '#F4A460', '#FF3030',
    '#FFEBCD', '#B0C4DE', '#00CDCD', '#C0FF3E', '#FFD700', '#8B4513', '#4EEE94', '#CD3278',
    '#00E5EE', '#E3A869', '#CD853F', '#ADD8E6', '#CD2990', '#EEE5DE', '#66CD00', '#7B68EE',
    '#FFA54F', '#A2B5CD', '#BC8F8F', '#8B2323', '#EE30A7', '#EEEED1', '#AEEEEE', '#5E2612',
    '#FF7F00', '#FFC0CB', '#EE3B3B', '#9370DB', '#848484', '#292421', '#CDBA96', '#B4EEB4',
    '#40E0D0', '#8B795E', '#3D9140', '#CDB7B5', '#CAE1FF', '#F0FFFF', '#2E8B57', '#FF6103',
    '#87CEEB', '#CD00CD', '#CDAA7D', '#836FFF', '#EEB4B4', '#8B7355', '#F0E68C', '#CDCDB4',
    '#B4CDCD', '#F0FFF0', '#00EEEE', '#708090', '#9AFF9A', '#FFA07A', '#FFB5C5', '#00688B',
    '#8A3324', '#191970', '#308014', '#FF83FA', '#838B8B', '#808A87', '#00FF7F', '#FFA500',
    '#EEAD0E', '#CD3333', '#4876FF', '#7CCD7C', '#EE5C42', '#AAAAAA', '#DAA520', '#8B3A3A',
    '#FFFAF0', '#B2DFEE', '#00EE76', '#FFFAFA', '#800080', '#C5C1AA', '#EEE685', '#FF3E96',
    '#EE0000', '#FDF5E6', '#EECFA1', '#8DB6CD', '#FF7256', '#7CFC00', '#838B83', '#BF3EFF',
    '#8B6914', '#00CD66', '#A4D3EE', '#00868B', '#8DEEEE', '#8B1C62', '#CDBE70', '#9F79EE', 
    '#C1CDC1', '#CD69C9', '#E0EEEE', '#8B7E66', '#8A2BE2', '#CDCD00', '#97FFFF', '#EEAEEE', 
    '#DC143C', '#CD919E', '#528B8B', '#CD6889', '#E6E6FA', '#E3CF57', '#4B0082', '#FF9912',
    '#F0F8FF', '#FF7F50', '#6CA6CD', '#8B8B83', '#F4F4F4', '#548B54', '#48D1CC', '#C1CDCD', 
    '#E0EEE0', '#3D59AB', '#FFB90F', '#FFD39B', '#8B5A2B', '#9C661F', '#EEE9BF', '#BCEE68',
    '#8EE5EE', '#8B0A50', '#FFF68F', '#EEA2AD', '#CD5B45', '#7FFF00', '#8B8378', '#9BCD9B',
    '#EEE8AA', '#8E8E38', '#668B8B', '#B3EE3A', '#00C78C', '#FFC125', '#8B475D', '#D8BFD8',
    '#FFE4C4', '#96CDCD', '#CDB5CD', '#00C5CD', '#00CED1', '#008B00', '#B8860B', '#1C86EE',
    '#EEC591', '#E066FF', '#B7B7B7', '#DEB887', '#FF6EB4', '#6959CD', '#90EE90', '#8B4789',
    '#EE7AE9', '#8968CD', '#D2B48C', '#FFFFE0', '#CDC9C9', '#BDB76B', '#00C957', '#EEDC82',
    '#3CB371', '#F5FFFA', '#B9D3EE', '#F5F5DC', '#0000CD', '#FF8247', '#EED5B7', '#FFEC8B',
    '#EE7600', '#7171C6', '#8B636C', '#8B814C', '#FFE4B5', '#1E1E1E', '#4F94CD', '#CDAD00', 
    '#CD5555', '#71C671', '#8B7500', '#473C8B', '#B0E0E6', '#FFFF00', '#8B4C39', '#006400',
    '#53868B', '#8B2252', '#FFB6C1', '#63B8FF', '#FFAEB9', '#EE6A50', '#87CEFF', '#87CEFA',
    '#5B5B5B', '#ADFF2F', '#008B45', '#EE4000', '#8A360F', '#8B6969', '#00008B', '#DB7093',
    '#7EC0EE', '#EE799F', '#CD6090', '#C76114', '#8B8682', '#458B74', '#FFF5EE', '#76EE00',
    '#000080', '#228B22', '#8B8B00', '#CD950C', '#EE82EE', '#282828', '#F5DEB3', '#3A5FCD',
    '#00FA9A', '#C67171', '#D1EEEE', '#8B5742', '#8B3E2F', '#CD3700', '#9AC0CD', '#555555',
    '#8B8989', '#EED8AE', '#551A8B', '#778899', '#FFFACD', '#458B00', '#008000', '#FFFFF0',
    '#EEB422', '#5CACEE', '#CD4F39', '#CDC0B0', '#FF7D40', '#8E388E', '#6E7B8B', '#CDC673',
    '#7A378B', '#E0FFFF', '#FFFFFF', '#6C7B8B', '#FFC1C1', '#8B4726', '#515151', '#CD9B1D',
    '#FF6347', '#FF34B3', '#FF0000', '#B0E2FF', '#8B3A62', '#CD5C5C', '#A2CD5A', '#00EE00',
    '#FF6A6A', '#CD6600', '#FFEFDB', '#E9967A', '#EEE9E9', '#7A67EE', '#CD8162', '#00F5FF',
    '#FFEFD5', '#CDAF95', '#00BFFF', '#CDB79E', '#1E90FF', '#EE2C2C', '#8B6508', '#FF7F24',
    '#8FBC8F', '#66CDAA', '#6495ED', '#EEE0E5', '#C1C1C1', '#B22222', '#EE00EE', '#FF82AB',
    '#AB82FF', '#79CDCD', '#7D26CD', '#03A89E', '#8B008B', '#5D478B', '#8B3626', '#808069',
    '#FFE4E1', '#EEDFCC', '#9400D3', '#BFEFFF', '#8B7D6B', '#FF8C69', '#C6E2FF', '#FF4500',
    '#FFE7BA', '#872657', '#808080', '#EE9572', '#CD8500', '#8B5A00', '#9932CC', '#EECBAD',
    '#CD8C95', '#CD1076', '#7D9EC0', '#104E8B', '#8B668B', '#698B22', '#EEE8CD', '#DDA0DD',
    '#4169E1', '#DA70D6', '#DCDCDC', '#68228B', '#CDC8B1', '#000000', '#6B8E23', '#FF69B4',
    '#800000', '#5F9EA0', '#8B4500', '#FCE6C9', '#D3D3D3', '#CDB38B', '#607B8B', '#F08080',
    '#CD9B9B', '#76EEC6', '#FAEBD7', '#68838B', '#EAEAEA', '#7FFFD4', '#C0C0C0', '#EE9A49',
    '#4A708B', '#008080', '#7AC5CD', '#98F5FF', '#8B2500', '#FFF0F5', '#8B8970', '#8B8878',
    '#6A5ACD', '#4682B4', '#EEEEE0', '#27408B', '#00FF00', '#FFDEAD', '#CD2626', '#CD96CD',
    '#9B30FF', '#36648B', '#F8F8FF', '#EEC900', '#EEEE00', '#FFE1FF', '#C1FFC1', '#CDC5BF',
    '#A0522D', '#8B5F65', '#CDC1C5', '#EE7621', '#FFBBFF', '#CD6839', '#698B69', '#BDFCC9',
    '#CD661D', '#FAFAD2', '#CDCDC1', '#FFF8DC', '#B452CD', '#8E8E8E', '#8470FF', '#483D8B',
    '#BBFFFF', '#0000FF', '#EE6AA7', '#EE7942', '#00CD00', '#9ACD32', '#C71585', '#EE9A00',
    '#CAFF70', '#32CD32', '#8B0000', '#B0171F', '#98FB98', '#8B1A1A', '#00B2EE', '#20B2AA',
    '#009ACD', '#A52A2A', '#EE6363', '#FAF0E6', '#8B7D7B', '#9A32CD', '#FF8000', '#7A8B8B',
    '#CD7054', '#9FB6CD', '#CDC9A5', '#D02090', '#00FFFF', '#CD0000', '#43CD80', '#FA8072',
    '#FFDAB9', '#D2691E']
shuffle(colors)
shuffle(colors)
# ########################################################################################################
# Some python functions###################################################################################
# ########################################################################################################
def load_cva_LDs(prefix):
	''' Read and load the CVA file '''
	f = open(prefix+'.cva')
	d={}
	for line in f:
		if line.startswith('""'):
			continue
		else:
			bl = line.strip().split(';')
			k=bl[0].split('"')[1]
			if k not in d:
				d[k]=[[],[]]#,[]]
			for i in range(2):
				d[k][i].append(bl[i+1].split('"')[1])
	return d

def plot2d(prefix,d,fontsize,symlog):
	markers = ['k.','b+','g*','r.','c+','m*','y.','k+','b*','g.','r+','c*','m.','y+','k*','b.',
	           'g+','r*','c.','m+','y*']
	ann=[]
	fig = plt.figure()
	ax = fig.add_subplot(111)    
	ax.spines['top'].set_color('none')
	ax.xaxis.tick_bottom()
	ax.spines['right'].set_color('none')
	ax.yaxis.tick_left()
	co=0
	ck=0
	for k in d:
		ck+=1
		x = d[k][0]
		y = d[k][1]
		#z = d[k][0][2]
		try:
			g=int(k)
			xt=map(float,x)
			yt=map(float,y)
			#nx, min_max_x, cx, vx, skewx, kurtx = stats.describe(xt)
			#ny, min_max_y, cy, vy, skewy, kurty = stats.describe(yt)
			cx=mean(xt)
			cy=mean(yt)
			ann.append((k,(float(cx)+0.1,float(cy)+0.1),colors[co]))
			if len(x) == 1:
				vx=0
				vy=0
				width= 0
				height= 0
			else:
				vx=var(xt)
				vy=var(yt)
				#ciy = stats.norm.interval(0.05,loc=cy,scale=sqrt(vy))
				#cix = stats.norm.interval(0.05,loc=cx,scale=sqrt(vx))
				width= 2*vx#cix[1] - cix[0]
				height=2*vy#ciy[1] - ciy[0]
			c = Ellipse((cx,cy),width,height)
							#radius=ci,#(max(max(xt),max(yt))-min(min(xt),min(yt)))+max(vx,vy),	                fill=False, label=k )
			ax.add_artist(c)
			c.set_clip_box(ax.bbox)
			c.set_edgecolor( colors[co] )
			c.set_facecolor( 'none' )  # "none" not None
			c.set_alpha( 1 )
			ax.plot(x,y, marker='o',  markersize=6, label=k, 
			        lw = 0, color=colors[co])
			ax.annotate(k,(float(cx)+(2*vx),float(cy)+2*(vy)),color=colors[co],fontsize=fontsize)
			co+=1
		except:
			if len(d.keys())-1 == ck:
				ax.plot(x,y, marker='o',  label=k,#markersize=6,
				        lw = 0 , color=colors[co])#color='0.95',label='Not classified')#disabled temporarily for gp120... get back if cazy
			else:
				ax.plot(x,y, marker='o',  #markersize=6, 
				        label=k, lw = 0 ,color=colors[co])#color='0.85')
				handles, labels = ax.get_legend_handles_labels()
				ax.legend(handles[::-1], labels[::-1])
			co+=1
	if symlog:
		ax.set_xscale("symlog")
		ax.set_yscale("symlog")
		ax.set_xlabel('LD 1 (symmetrical log)', fontsize=fontsize)
		ax.set_ylabel('LD 2 (symmetrical log)', fontsize=fontsize)
	else:
		ax.set_xlabel('LD 1', fontsize=fontsize)
		ax.set_ylabel('LD 2', fontsize=fontsize)
	#ax.set_aspect('equal')
	'''ax.legend(#bbox_to_anchor=(0.45, -0.2), 
              ncol=5+(co/5),
              loc=0, 
              prop={'size':10},fancybox=True, shadow=True)'''
	fig.tight_layout()
	plt.show()
	fig.savefig(prefix+'_CVA2D.png', dpi=300)

# ######################################################################################################
# Aplication of the code ###############################################################################
# ######################################################################################################
if __name__ == "__main__":
	if len(sys.argv) == 1 or '-help' in sys.argv:
		print 'Usage: python PlotCVA.py <prefix> [options]'
		print 'A <prefix>.cva must be provided, and comming from R. First row are headers.'\
		      ' First column Group of the CVA, second first LD, third second LD, and fourth'\
		      ' the third LD. All entries are between double quotes (").'
		print '[Options]:'
		print '-3d : Use this if you want to plot the 3 dimensions. ( Default : No )'
		print '\t-symlog : Use symmetrical log transformation of the axes. In only positive'\
		      ' values behaves like log. ( Default: No )'
		print '\t-fontsize=XX : Will use XX as the font size for axes titles and numbering'\
		      ' in the plot. ( Default : 12 )'

	# Default Parameters ###############################################################################
	prefix = sys.argv[1]
	symlog=False
	threeD = False
	# Get user input ###################################################################################
	for arg in sys.argv[1:]:
		if arg == '-symlog':
			symlog=True
		elif arg.startswith('-fontsize='):
			fontsize= int(arg[10:])
		elif arg == '-3D':
			threeD = True


	d=load_cva_LDs(prefix)
	#if threeD:
	plot2d(prefix,d,fontsize,symlog)