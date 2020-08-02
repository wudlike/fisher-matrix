import optimization_value as opt
#import test3 as opt
import time
import argparse
import settings as sts 
import numpy as np

parser = argparse.ArgumentParser(description='test running')
parser.add_argument('-Ap','--Add_prior',default='True', choices=['True','False'],
                    help='Whether to consider prior information, default=True')
parser.add_argument('-fxn','--Fixed_Ndet',default='False', choices=['True','False'],
                    help='Fixed detector numbers, default=False')
parser.add_argument('-skyf','--sky_fraction',default='large',choices=['small','large'],
                    help='Choose a small(0.05)/large(0.7) sky fraction, default=small')
parser.add_argument('-r','--r',default=0.01,choices=['0','0.1','0.01'],
                    help='input tensor to scalar ratio, default=0.01')
parser.add_argument('-nu1',default=95,type=int,
                    help='input first frequency, default=95')
parser.add_argument('-nu2',default=False,type=int,
                    help='input second frequency, default=False')
parser.add_argument('-nu3',default=False,type=int,
                    help='input third frequency, default=False')
parser.add_argument('-nu4',default=False,type=int,
                    help='input fourth frequency, default=False')
parser.add_argument('-nu5',default=False,type=int,
                    help='input fourth frequency, default=False')
args = parser.parse_args()

config = vars(args)
class Switch():
    def __init__(self):
        self.skyf=args.sky_fraction
        self.r = args.r
        self.nu1 = args.nu1
        self.nu2 = args.nu2
        self.config = vars(args)
        
    def test(self):
        return 0

if __name__ == '__main__':
    st = sts.Settings()
    start = time.time()
    sts.test(config)
    print("running...")
    start = time.time()
    if config['Fixed_Ndet']=='True':
        if args.nu3 == 0:
            num_x,num_y,sigma_min_rr=opt.Fisher(args.Add_prior,config).optimize(config)
        elif args.nu4 == 0:
            num_x,num_y,num_z,sigma_min_rr=opt.Fisher(args.Add_prior,config).optimize(config)    
        elif args.nu5 == 0:
            num_x,num_y,num_z,num_m,sigma_min_rr=opt.Fisher(args.Add_prior,config).optimize(config)  
        elif args.nu5 != 0:
            num_x,num_y,num_z,num_m,num_n,sigma_min_rr=opt.Fisher(args.Add_prior,config).optimize(config)
    else:
        if st.n_nu == 4:
            sigma_r_list,num_x,num_y,num_z,num_m,x_list,y_list,z_list,m_list,sigma_rr_list, sigma_min_rr, sigmaMat_min,sigmaMat_sqrt_diag = opt.Fisher(args.Add_prior,config).optimize(config)    
        else:
            sigma_r_list,num_x,num_y,num_z,x_list,y_list,z_list,sigma_rr_list, sigma_min_rr, sigmaMat_min,sigmaMat_sqrt_diag = opt.Fisher(args.Add_prior,config).optimize(config)

    # save files
#    opt.Fisher(args.Add_prior).save_file(config)            
    end = time.time()
    print('\n''################ time cost ###############')
    running_time = (end - start)/60    #minutes
    print('running time: %.5f mins' %running_time) 
    if config['Fixed_Ndet']=='True':
        print('sigma_min_rr= %.9f '%np.array(sigma_min_rr[0]))
        print('################## answer ################')
        print('sigma_min_rr= %.5f '%sigma_min_rr[0])
    else:
        print('sigma_min_rr= %.9f '%np.array(sigma_min_rr))
        print('################## answer ################')
        print('sigma_min_rr= %.5f '%sigma_min_rr)
    if args.nu2 == 0:
        print(str(args.nu1)+' GHz = %d '%num_x)
    if args.nu2 != 0 and args.nu3 == 0:
        print(str(args.nu1)+' GHz = %d '%num_x)
        print(str(args.nu2)+' GHz = %d '%num_y)
    if args.nu3 != 0 and args.nu4 ==0:
        print(str(args.nu1)+' GHz = %d '%num_x)
        print(str(args.nu2)+' GHz = %d '%num_y)
        print(str(args.nu3)+' GHz = %d '%num_z)
    if args.nu3 != 0 and args.nu4 != 0 and args.nu5 == 0:
        print(str(args.nu1)+' GHz = %d '%num_x)
        print(str(args.nu2)+' GHz = %d '%num_y)
        print(str(args.nu3)+' GHz = %d '%num_z)
        print(str(args.nu4)+' GHz = %d '%num_m)
    if args.nu5 != 0:
        print(str(args.nu1)+' GHz = %d '%num_x)
        print(str(args.nu2)+' GHz = %d '%num_y)
        print(str(args.nu3)+' GHz = %d '%num_z)
        print(str(args.nu4)+' GHz = %d '%num_m)
        print(str(args.nu5)+' GHz = %d '%num_n)        

        
        
        
        
        
        
        
