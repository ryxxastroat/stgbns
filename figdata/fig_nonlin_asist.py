import sys
import timeit
import shutil
import os
import numpy as np

t0 = timeit.time.time()

def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

def piecewisef(x):
   return np.piecewise(x, [1<= x <=5, 6<= x <=10, 11<= x <=15, 16<= x <=20, 21<= x <=25, 26<= x <=30, 31<= x <=35, 36<= x <=40], [5, 6, 7, 8, 1, 2, 3, 4])


# setup the parameters
parameter1 = ['0.1**2', '1.**2']
parameter2 = ['-1.', '-0.2', '3.', '10.']
parameter3 = ['WFF1','SLy4','AP4','MPA1','PAL1']
ip, jp, kp = len(parameter1), len(parameter2), len(parameter3)
 
intset = np.arange(1, ip*jp*kp+1, dtype=int)
#print(intset)


# create files 
f = open("test.py", 'w+')


for i in range(0, ip):
   for j in range(0, jp):
      for k in range(0, kp):
         dataname = 'stdata' +  str(i*jp*kp+j*kp+k+1)         
         if i*jp*kp+j*kp+k+1 !=29:
           f.write(dataname + '= np.genfromtxt(\'stgb_tid_v1_comb_data' + str(i*jp*kp+j*kp+k+1) + '.txt\') \n')
           f.write('c0'+str(i*jp*kp+j*kp+k+1)+', c1'+str(i*jp*kp+j*kp+k+1)+', c2'+str(i*jp*kp+j*kp+k+1)+', c3'+str(i*jp*kp+j*kp+k+1)+', c4'+str(i*jp*kp+j*kp+k+1)+', c5'+str(i*jp*kp+j*kp+k+1)+', c6'+str(i*jp*kp+j*kp+k+1)+', c7'+str(i*jp*kp+j*kp+k+1)+', c8'+str(i*jp*kp+j*kp+k+1)+', c9'+str(i*jp*kp+j*kp+k+1)+', c10'+str(i*jp*kp+j*kp+k+1)+'='+dataname+'[:, 0], '+dataname+'[:, 1], '+dataname+'[:, 2], '+dataname+'[:, 3], '+dataname+'[:, 4], '+dataname+'[:, 5], '+dataname+'[:, 6], '+dataname+'[:, 7], '+dataname+'[:, 8], '+dataname+'[:, 9], '+dataname+'[:, 10] \n')
           f.write('x'+str(i*jp*kp+j*kp+k+1)+' = c3'+str(i*jp*kp+j*kp+k+1)+'/MSUN \n')  
           f.write('y'+str(i*jp*kp+j*kp+k+1)+' = c4'+str(i*jp*kp+j*kp+k+1)+' \n') 
            
           f.write('ax'+str(piecewisef(i*jp*kp+j*kp+k+1))+'.plot(x'+str(i*jp*kp+j*kp+k+1)+', y'+str(i*jp*kp+j*kp+k+1)+', \'o\', color = colorset['+str( (i*jp*kp+j*kp+k)%5 )+']) \n')  
           #f.write('ax'+str(piecewisef(i*jp*kp+j*kp+k+1))+'.semilogy(x'+str(i*jp*kp+j*kp+k+1)+', y'+str(i*jp*kp+j*kp+k+1)+', color = colorset['+str( (i*jp*kp+j*kp+k)%5 )+']) \n')                   

           
f.close  



print( '\n *** STG_solver uses %.2f seconds\n' % (timeit.time.time() - t0))
