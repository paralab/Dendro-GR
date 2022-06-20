import numpy as np

# total mass
M=1
N=40
R=800
# mass ratio
q=np.array([1       ,2      ,4      ,8      ,16     ,32     ,64     ,128])


#number of elements, 
gEsz=np.array([ 21624, #1
               37696, #2
               60096, #4
               23976, #8
               39152, #16
               39936, #32
               194468,#64
               1321356])#128


#cg nodes
gNsz=np.array([24936973, #1
               34731915, #2 
               48071805, #4 
               26243162, #8
               36153655, #16
               35880015, #32
               109821531, #64
               495862508]) #128

m1=M/(q+1)
m2=M-m1

# radius approx. SR
r1=0.5*m1
r2=0.5*m2

x1=2*r1/N
x2=2*r2/N

d1=np.floor(np.log2((R/(6*x1))) + 2)
d2=np.floor(np.log2((R/(6*x2))) + 2)
#print(m1*2.0)
#print(m2*2.0)
#print(x2/x1)
amr1=m1*2
amr2=m2*2
cw=15
import texttable as tt
tab = tt.Texttable()
headings =         ['q','m1','m2','dx1','dx2','lev 1', 'lev 2','amr r1','amr r2','octants','grid pts']
tab.set_cols_dtype(['i','e' ,'e' ,'e'  ,'e'  ,'i'    , 'i'    ,'e'     ,'e'     , 'i'     ,'e'])
tab.set_cols_width([cw ,cw  ,cw  , cw  , cw  ,cw     , cw     ,cw      ,cw      , cw      , cw ])
tab.header(headings)


for r in range(0,len(q)):
	row=[q[r],m1[r],m2[r],x1[r],x2[r],d1[r],d2[r],amr1[r],amr2[r],gEsz[r],gNsz[r]]
	#print(row)
	tab.add_row(row)
	
s = tab.draw()
print (s)
