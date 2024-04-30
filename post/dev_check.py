import numpy as np
import sys
import time

## Settings
directory1 = '/home/mani/kimliu/slip-b151/output/flow_field/'
#directory2 = '/home/mani/kimliu/pat-L077/output/flow_field_transportee/'
#directory3 = '/home/mani/kimliu/mfm-01/output/flow_field_IMFM/'

Nx, Ny, Nz = 192, 128, 192
Lx, Ly, Lz = 2.*np.pi, 2., 1.*np.pi
Np = Nx*Ny*Nz
Nh = int(Ny/2)
Re = 197.5


startT=float(sys.argv[1])
endT=float(sys.argv[2])
dT=float(sys.argv[3])

start = time.time()
print('Time: ' + repr(startT) + ' ~ ' + repr(endT) + ', ' + repr(dT))


## Setting t
Nt=round((endT-startT)/dT)+1
t_filename=np.linspace(startT,endT,Nt)
print(t_filename)

## Setting coordinates
xc=np.linspace(Lx/(2*Nx), Lx-Lx/(2*Nx), Nx, endpoint=True)
xf = xc - Lx/(2*Nx)
zc=np.linspace(Lz/(2*Nz), Lz-Lz/(2*Nz), Nz, endpoint=True)
zf = zc - Lz/(2*Nz)

gamma=2.
jf=np.linspace(1, Ny-1, Ny-1, endpoint=True)
yf=-np.tanh(gamma*(1.-2./Ny*jf))/np.tanh(gamma)
jc=np.linspace(1, Ny, Ny, endpoint=True)
yc=-np.tanh(gamma*(1.-2./Ny*(jc-0.5)))/np.tanh(gamma)
yf = np.append(yf,np.array([1])) #face value (Ny-1) + endpoint located on top wall


## Initialize arrays
t = np.ndarray(shape=(Nt))
ucl = np.ndarray(shape=(Nt))
ub = np.ndarray(shape=(Nt))
ut = np.ndarray(shape=(Nt))
u1_prof = np.ndarray(shape=(Ny,Nt))


# e.g. field_name(30.2) returns "field_030_200.dat"
def field_name(t):
    print(t)
    i = int(t)
    if i == 1000:
        return "field_***_000.dat"
    f = round(t - i, 3)
    while f >= 1:
      i = i + 1
      f = round(t-i,3)
    if f != 0 :
        return "field_%s_%s.dat" % (str(i).zfill(4), str(int(f*1000)).zfill(3))
    if f == 0 :
        return "field_%s_%s.dat" % (str(i).zfill(4), '000')


## Iterate through files
for f in range(Nt):
    # Load u flow field
    loc = 0
    loc += np.int32().itemsize # io_ver

    fo = open(directory1+field_name(t_filename[f]),'rb')
    fo.seek(loc)
    a_f = np.fromfile(fo,dtype=np.float64,count=1)
    t[f] = np.array(a_f[0])

    loc += np.float64().itemsize
    fo.seek(loc)
    a_f = np.fromfile(fo,dtype=np.float64,count=Np)
    u1 = np.array(a_f[0:Np])

    fo.close()
 
    u1 = np.reshape(u1, (Nx, Ny, Nz), order='F')
    u1m = np.mean(u1,axis=(0,2))
     
    # Post process

    # calculate us
    ub[f] = u1m[0]
    ut[f] = u1m[-1]
    ucl[f] = (u1m[Nh]+u1m[Nh-1])/2.0
    u1_prof[:,f] = u1m
     
# Save
np.savez('data/donor_dev', tdev=t, ucl=ucl, ub=ub, ut=ut, u1m=u1_prof)

print('Success!')
print(time.time()-start)
