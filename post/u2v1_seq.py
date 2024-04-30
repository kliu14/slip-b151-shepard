import numpy as np
import sys
import time

## Settings
directory1 = '/home/mani/kimliu/slip-shep/output/flow_field/'
directory2 = '/home/mani/kimliu/slip-shep/output/flow_field_transportee/'
#directory3 = '/home/mani/kimliu/mfm-01/output/flow_field_IMFM/'

Nx, Ny, Nz = 192, 128, 192
Lx, Ly, Lz = 2.*np.pi, 2., 1.*np.pi
Np = Nx*Ny*Nz
Re = 197.5

startT=float(sys.argv[1])
endT=float(sys.argv[2])
dT=float(sys.argv[3])
nn = int(endT)-int(startT)
Nl = int(np.log2(nn)) + 1

start = time.time()
print('Time: ' + repr(startT) + ' ~ ' + repr(endT) + ', ' + repr(dT))


## Setting t
Nt=round((endT-startT)/dT)+1
#t_filename=np.linspace(startT,endT,Nt)
t_filename = np.linspace(endT,startT,Nt)
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

# variables for u
u1_txz = np.ndarray(shape=(Ny))
u2_txz = np.ndarray(shape=(Ny))
u3_txz = np.ndarray(shape=(Ny))

u1_txz.fill(0.)
u2_txz.fill(0.)
u3_txz.fill(0.)

# variables for v
v1_txz = np.ndarray(shape=(Ny))
v2_txz = np.ndarray(shape=(Ny))
v3_txz = np.ndarray(shape=(Ny))

v1_txz.fill(0.)
v2_txz.fill(0.)
v3_txz.fill(0.)

# uv's located in y_f
u2v1 = np.ndarray(shape=(Nx,Ny-1,Nz))
u2v1_txz = np.ndarray(shape=(Ny-1))
u2u1 = np.ndarray(shape=(Nx,Ny-1,Nz))
u2u1_txz = np.ndarray(shape=(Ny-1))

u2v1.fill(0.)
u2v1_txz.fill(0.)
u2u1.fill(0.)
u2u1_txz.fill(0.)


# e.g. field_name(30.2) returns "field_030_200.dat"
def field_name(t):
    i = int(t)
    if i == 1000:
        return "field_***_000.dat"
    f = round(t - i, 3)
    if f != 0 :
        return "field_%s_%d.dat" % (str(i).zfill(4), int(f*1000))
    if f == 0 :
        return "field_%s_%s.dat" % (str(i).zfill(4), '000')


## Iterate through files

for g in range(Nl):
    Nt = int((2**g)/dT)
    Nt_prev = 0

    if g != 0:
        Nt_prev = int(Nt/2)
        u1_txz = u1_txz / 2.0
        u2_txz = u2_txz / 2.0
        u3_txz = u3_txz / 2.0
        v1_txz = v1_txz / 2.0
        v2_txz = v2_txz / 2.0
        v3_txz = v3_txz / 2.0
        u2v1_txz = u2v1_txz / 2.0
        u2u1_txz = u2u1_txz / 2.0

    for f in range(Nt_prev,Nt):
        # Load u flow field
        loc = 0
        loc += np.int32().itemsize # io_ver

        fo = open(directory1+field_name(t_filename[f]),'rb')
        fo.seek(loc)
        #a_f = np.fromfile(fo,dtype=np.float64,count=1)
        #t[f] = np.array(a_f[0])

        loc += np.float64().itemsize
        fo.seek(loc)
        a_f = np.fromfile(fo,dtype=np.float64,count=Np)
        u1 = np.array(a_f[0:Np])

        loc += Np*np.float64().itemsize
        fo.seek(loc)
        a_f = np.fromfile(fo,dtype=np.float64,count=Np)
        u2 = np.array(a_f[0:Np])

        loc += Np*np.float64().itemsize
        fo.seek(loc)
        a_f = np.fromfile(fo,dtype=np.float64,count=Np)
        u3 = np.array(a_f[0:Np])

        fo.close()
     
        u1 = np.reshape(u1, (Nx, Ny, Nz), order='F')
        u2 = np.reshape(u2, (Nx, Ny, Nz), order='F')
        u3 = np.reshape(u3, (Nx, Ny, Nz), order='F')

        # Load v flow field
        loc = 0
        loc += np.int32().itemsize # io_ver

        fo = open(directory2+field_name(t_filename[f]),'rb')
        fo.seek(loc)
        #a_f = np.fromfile(fo,dtype=np.float64,count=1)
        #t[f] = np.array(a_f[0])

        loc += np.float64().itemsize
        fo.seek(loc)
        a_f = np.fromfile(fo,dtype=np.float64,count=Np)
        v1 = np.array(a_f[0:Np])

        loc += Np*np.float64().itemsize
        fo.seek(loc)
        a_f = np.fromfile(fo,dtype=np.float64,count=Np)
        v2 = np.array(a_f[0:Np])

        loc += Np*np.float64().itemsize
        fo.seek(loc)
        a_f = np.fromfile(fo,dtype=np.float64,count=Np)
        v3 = np.array(a_f[0:Np])

        fo.close()
     
        v1 = np.reshape(v1, (Nx, Ny, Nz), order='F')
        v2 = np.reshape(v2, (Nx, Ny, Nz), order='F')
        v3 = np.reshape(v3, (Nx, Ny, Nz), order='F')


        # Post process

        #calculate u2v1 and u2u1
        for i in range(Nx):
            for j in range(Ny-1):
                u2v1[i,j,:] = (u2[i-1,j,:]+u2[i,j,:])*(v1[i,j,:]+v1[i,j+1,:])/4.0
                u2u1[i,j,:] = (u2[i-1,j,:]+u2[i,j,:])*(u1[i,j,:]+u1[i,j+1,:])/4.0
                 
        # add mean u
        u1_txz = u1_txz + np.mean(u1, axis=(0, 2)) / Nt
        u2_txz = u2_txz + np.mean(u2, axis=(0, 2)) / Nt
        u3_txz = u3_txz + np.mean(u3, axis=(0, 2)) / Nt

        # add mean v
        v1_txz = v1_txz + np.mean(v1, axis=(0, 2)) / Nt
        v2_txz = v2_txz + np.mean(v2, axis=(0, 2)) / Nt
        v3_txz = v3_txz + np.mean(v3, axis=(0, 2)) / Nt

        # add mean u2v1 and u2u1
        u2v1_txz = u2v1_txz + np.mean(u2v1,axis=(0,2)) / Nt
        u2u1_txz = u2u1_txz + np.mean(u2u1,axis=(0,2)) / Nt

        print(t_filename[f])
   
    # save
    np.savez('data/u2v1_t{:03d}'.format(int(Nt*0.2)), yc=yc, yf=yf[:-1], u1=u1_txz, u2=u2_txz[:-1], u3=u3_txz, v1=v1_txz, v2=v2_txz[:-1], v3=v3_txz, uv=u2v1_txz, uu=u2u1_txz) 
    print('Averaging window ' + repr(Nt*0.2) + ' done!') 
    print(time.time()-start)

print('Success!')
print(time.time()-start)
