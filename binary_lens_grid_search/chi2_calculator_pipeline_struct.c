//#include<iostream>
//#include<fstream>
#include<math.h>
//#include<stdbool.h>
//#include<string>
//#include<Python.h>
//#include<time.h>
//#include<stdio.h>

/*float checking_whether_densed_to_next_layer(float x,float y,bool arg_x,bool arg_y,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,int sequence_number_in_next_layer_file,int current_layer)
{
	float A ;

	x = 2 * x - arg_x ;
	y = 2 * y - arg_y ;
	
	int arg_x_plus_double_arg_y = 2 * arg_y + arg_x ;

	int arg = layer_length_accumulated_in_front_of_each_layer[current_layer] + 4 * sequence_number_in_next_layer_file + arg_x_plus_double_arg_y ;
	
	if ( all_layer_whether_densed[arg] )
	{
		arg_x = int( 2 * x ) ;
		arg_y = int( 2 * y ) ;

		sequence_number_in_next_layer_file = all_layer_sequence_number_in_next_layer_file[arg] ;
        // must wait for 50 ns here
		A = checking_whether_densed_to_next_layer(x,y,arg_x,arg_y,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,sequence_number_in_next_layer_file,current_layer+1) ;
	}
	else
	{
		float A1 , A2 , A3 , A4 ;
		A1 = all_layer_corner_mag[ 4 * arg ] ;
		A2 = all_layer_corner_mag[ 4 * arg + 1 ] ;
		A3 = all_layer_corner_mag[ 4 * arg + 2 ] ;
		A4 = all_layer_corner_mag[ 4 * arg + 3 ] ;
		
		A = (1-x)*(1-y)*A1 + x*(1-y)*A2 + (1-x)*y*A3 + x*y*A4 ;
	}

	return A ;
}*/

/*void interpolating_to_get_N_magnification(float *x_list,float *y_list,float *A_list,int Npoint,\
                                        register bool Whether0,register int Offset0,\
                                        register bool Whether1,register int Offset1,\
                                        register bool Whether2,register int Offset2,\
                                        register bool Whether3,register int Offset3,\
                                        register bool Whether4,register int Offset4,\
                                        register bool Whether5,register int Offset5,\
                                        register bool Whether6,register int Offset6,\
                                        register bool Whether7,register int Offset7,\
                                        float *all_layer_corner_mag,\
                                        bool *all_layer_whether_densed,\
                                        int *all_layer_sequence_number_in_next_layer_file,\
                                        int *layer_length_accumulated_in_front_of_each_layer,\
                                        float box_size)
{
    // these N points must be in the map
    for(int i=0;i<Npoint;i++)
        {
            Whether[i] = all_layer_whether_densed[arg] ;
            Offset[i]  = all_layer_sequence_number_in_next_layer_file[arg] ;
        }
}*/




struct whether_offset_4mag
{
        _Bool whether_densed ;
        int sequence_number_in_next_layer_file ;
        float corner_mag_1 ;
        float corner_mag_2 ;
        float corner_mag_3 ;
        float corner_mag_4 ;
};

// fix Npoint=8 then compile, then the compiler can do more optimization
// the resulting speed is 1.65 s compared to 1.9 s of the flexiable Npoint version

// But with -march=skylake-avx512 and -I/usr/include/python2.7/ 
// the Npoint fixed to 8 is 1.64 s compared to 2.1 s of the flexiable Npoint version

///////// test in the environment that I occupy the whole node (192 thread) but run the code only on one thread /////////
///////// the result is then more stable                                                                        /////////
void interpolating_to_get_8_magnification(float * __restrict x_list,float * __restrict y_list,float * __restrict A_list,\
                                        struct whether_offset_4mag * __restrict one_map_all_layer_struct,\
                                        /*float *all_layer_corner_mag,\
                                        _Bool *all_layer_whether_densed,\
                                        int *all_layer_sequence_number_in_next_layer_file,\*/
                                        int * __restrict layer_length_accumulated_in_front_of_each_layer,\
                                        float box_size)
{
        // these 8 points must be in the map

        _Bool whether_list[8] ; 
        int offset_list[8] ; 

        int arg_list[8] ; 
        register int arg ;

        _Bool arg_x_list[8] ;
        _Bool arg_y_list[8] ;

        register int current_layer = 0 ;
        register int offset_in_front_of_current_layer ;

        _Bool finish_list[8] ; 

        for(register int i=0;i<8;i++) // vectorize
        {
                // set the initial value to be false
                finish_list[i] = 0 ;
        }

        int arg_retain_list[8] ; 
        float x_retain_list[8] ; 
        float y_retain_list[8] ; 
        register int n_already = 0 ;

        float A1_list[8] ; 
        float A2_list[8] ; 
        float A3_list[8] ; 
        float A4_list[8] ; 
        
        while ( n_already < 8 ) // if all A have been calculated, then jump out the while loop 
        {
                offset_in_front_of_current_layer = layer_length_accumulated_in_front_of_each_layer[current_layer] ;
                // this array will eventually be in the cache due to its small size

                if ( current_layer == 0 )
                {
                        for(register int i=0;i<8;i++) // vectorize
                        {
                                // the absolute coordinates in the map is adjusted to the relative coordinates in the layer 0's square
                                x_list[i] = 0.5 * ( x_list[i] + box_size ) / box_size ;
                                y_list[i] = 0.5 * ( y_list[i] + box_size ) / box_size ;

                                // calculate the layer 0's square's location in the DRAM
                                arg_list[i] = offset_in_front_of_current_layer + 4 * 0 + 0 ;
                        } 
                }
                else
                {
                        for(register int i=0;i<8;i++)
                        {       
                                // the new square is one of the four squares densed from the previous layer's big square
                                // calculate relative coordinates in the new square
                                x_list[i] = 2 * x_list[i] - arg_x_list[i] ;
                                y_list[i] = 2 * y_list[i] - arg_y_list[i] ;

                                // calculate the new square's location in the DRAM
                                // offset is the four squares' (densed from the previous layer's big square) location in the current layer
                                // (arg_x, arg_y) is the new square's index in these four squares
                                arg_list[i] = offset_in_front_of_current_layer + 4 * offset_list[i] + 2 * arg_y_list[i] +  arg_x_list[i] ;
                                // even if offset_list[i] is not in register, 
                                // it must be in the cache, because the array 'offset_list' will be called for many times
                        }
                }

                for(register int i=0;i<8;i++)
                {
                        // calculate the next-layer-square's index in next layer's four squares (assuming the new square is densed)
                        arg_x_list[i] = (int)( 2 * x_list[i] ) ;
                        arg_y_list[i] = (int)( 2 * y_list[i] ) ;
                }

                // keep the arg_list[i] unchanged for finish[i]==1 using arg_retain_list[i]
                // to insure the read operation in the RAM access 'for loop' is meaningful
                for(register int i=0;i<8;i++)
                {
                        if ( finish_list[i]==1 )
                        {
                                arg_list[i] = arg_retain_list[i] ;
                        }
                }

                for(register int i=0;i<8;i++) // hope cpu can do another memory access when waiting for one memory access
                // see https://stackoverflow.com/questions/45382914/parallel-memory-access-on-modern-processors
                {
                        arg = arg_list[i] ;

                        // 'whether' represents whether the new square is densed
                        // if the new square is densed :
                        //      'offset'  represents next layer's four squares' (densed from the new square) location in the next layer
                        // else :
                        //      'offset' = 65535 which is meaningless
                        whether_list[i] = one_map_all_layer_struct[arg].whether_densed ;
                        offset_list[i]  = one_map_all_layer_struct[arg].sequence_number_in_next_layer_file ;

                        // separate these two instructions from others
                        // hope to achieve instruction parallel and thus mask the latency of reading from DRAM

                        // Actually, no need to re-read the already finished cell's value
                        // Can use 'if ( finish_list[i]==0 )'
                        // But the resulting speed has no difference, because the already finished cell's value is already in the cache
                }

                for(register int i=0;i<8;i++)
                {
                        if ( whether_list[i]==0 && finish_list[i]==0 ) // to insure each point will only calculate magnification once
                        {
                                finish_list[i] = 1 ;

                                arg_retain_list[i] = arg_list[i] ;

                                x_retain_list[i] = x_list[i] ;
                                y_retain_list[i] = y_list[i] ;

                                n_already += 1 ;
                        }
                }

                current_layer += 1 ;

        }

        for(register int i=0;i<8;i++) // hope cpu can do another memory access when waiting for one memory access
        {
                arg = arg_retain_list[i] ;

                A1_list[i] = one_map_all_layer_struct[arg].corner_mag_1 ;
                A2_list[i] = one_map_all_layer_struct[arg].corner_mag_2 ;
                A3_list[i] = one_map_all_layer_struct[arg].corner_mag_3 ;
                A4_list[i] = one_map_all_layer_struct[arg].corner_mag_4 ;
                // these values are already in the cache because of the struct
        }

        for(register int i=0;i<8;i++) // vectorize
        {
                A_list[i] = (1-x_retain_list[i]) * (1-y_retain_list[i]) * A1_list[i] + \
                               x_retain_list[i]  * (1-y_retain_list[i]) * A2_list[i] + \
                            (1-x_retain_list[i]) *    y_retain_list[i]  * A3_list[i] + \
                               x_retain_list[i]  *    y_retain_list[i]  * A4_list[i] ;
        }

}




void interpolating_to_get_N_magnification(float *x_list,float *y_list,float *A_list,\
                                        register int Npoint,\
                                        struct whether_offset_4mag *one_map_all_layer_struct,\
                                        /*float *all_layer_corner_mag,\
                                        _Bool *all_layer_whether_densed,\
                                        int *all_layer_sequence_number_in_next_layer_file,\*/
                                        int *layer_length_accumulated_in_front_of_each_layer,\
                                        float box_size)
{
        // these N points must be in the map

        _Bool whether_list[Npoint] ; 
        int offset_list[Npoint] ; 

        int arg_list[Npoint] ; 
        register int arg ;

        _Bool arg_x_list[Npoint] ;
        _Bool arg_y_list[Npoint] ;

        register int current_layer = 0 ;
        register int offset_in_front_of_current_layer ;

        _Bool finish_list[Npoint] ; 

        for(register int i=0;i<Npoint;i++) // vectorize
        {
                // set the initial value to be false
                finish_list[i] = 0 ;
        }

        int arg_retain_list[Npoint] ; 
        float x_retain_list[Npoint] ; 
        float y_retain_list[Npoint] ; 
        register int n_already = 0 ;

        float A1_list[Npoint] ; 
        float A2_list[Npoint] ; 
        float A3_list[Npoint] ; 
        float A4_list[Npoint] ; 
        
        while ( n_already < Npoint ) // if all A have been calculated, then jump out the while loop 
        {
                offset_in_front_of_current_layer = layer_length_accumulated_in_front_of_each_layer[current_layer] ;
                // this array will eventually be in the cache due to its small size

                if ( current_layer == 0 )
                {
                        for(register int i=0;i<Npoint;i++) // vectorize
                        {
                                // the absolute coordinates in the map is adjusted to the relative coordinates in the layer 0's square
                                x_list[i] = 0.5 * ( x_list[i] + box_size ) / box_size ;
                                y_list[i] = 0.5 * ( y_list[i] + box_size ) / box_size ;

                                // calculate the layer 0's square's location in the DRAM
                                arg_list[i] = offset_in_front_of_current_layer + 4 * 0 + 0 ;
                        } 
                }
                else
                {
                        for(register int i=0;i<Npoint;i++)
                        {       
                                // the new square is one of the four squares densed from the previous layer's big square
                                // calculate relative coordinates in the new square
                                x_list[i] = 2 * x_list[i] - arg_x_list[i] ;
                                y_list[i] = 2 * y_list[i] - arg_y_list[i] ;

                                // calculate the new square's location in the DRAM
                                // offset is the four squares' (densed from the previous layer's big square) location in the current layer
                                // (arg_x, arg_y) is the new square's index in these four squares
                                arg_list[i] = offset_in_front_of_current_layer + 4 * offset_list[i] + 2 * arg_y_list[i] +  arg_x_list[i] ;
                                // even if offset_list[i] is not in register, 
                                // it must be in the cache, because the array 'offset_list' will be called for many times
                        }
                }

                for(register int i=0;i<Npoint;i++)
                {
                        // calculate the next-layer-square's index in next layer's four squares (assuming the new square is densed)
                        arg_x_list[i] = (int)( 2 * x_list[i] ) ;
                        arg_y_list[i] = (int)( 2 * y_list[i] ) ;
                }

                // keep the arg_list[i] unchanged for finish[i]==1 using arg_retain_list[i]
                // to insure the read operation in the RAM access 'for loop' is meaningful
                for(register int i=0;i<Npoint;i++)
                {
                        if ( finish_list[i]==1 )
                        {
                                arg_list[i] = arg_retain_list[i] ;
                        }
                }

                for(register int i=0;i<Npoint;i++) // hope cpu can do another memory access when waiting for one memory access
                // see https://stackoverflow.com/questions/45382914/parallel-memory-access-on-modern-processors
                {
                        arg = arg_list[i] ;

                        // 'whether' represents whether the new square is densed
                        // if the new square is densed :
                        //      'offset'  represents next layer's four squares' (densed from the new square) location in the next layer
                        // else :
                        //      'offset' = 65535 which is meaningless
                        whether_list[i] = one_map_all_layer_struct[arg].whether_densed ;
                        offset_list[i]  = one_map_all_layer_struct[arg].sequence_number_in_next_layer_file ;
                        // separate these two instructions from others
                        // hope to achieve instruction parallel and thus mask the latency of reading from DRAM
                }

                for(register int i=0;i<Npoint;i++)
                {
                        if ( whether_list[i]==0 && finish_list[i]==0 ) // to insure each point will only calculate magnification once
                        {
                                finish_list[i] = 1 ;

                                arg_retain_list[i] = arg_list[i] ;

                                x_retain_list[i] = x_list[i] ;
                                y_retain_list[i] = y_list[i] ;

                                n_already += 1 ;
                        }
                }

                current_layer += 1 ;

        }

        for(register int i=0;i<Npoint;i++) // hope cpu can do another memory access when waiting for one memory access
        {
                arg = arg_retain_list[i] ;

                A1_list[i] = one_map_all_layer_struct[arg].corner_mag_1 ;
                A2_list[i] = one_map_all_layer_struct[arg].corner_mag_2 ;
                A3_list[i] = one_map_all_layer_struct[arg].corner_mag_3 ;
                A4_list[i] = one_map_all_layer_struct[arg].corner_mag_4 ;
                // these values are already in the cache because of the struct
        }

        for(register int i=0;i<Npoint;i++) // vectorize
        {
                A_list[i] = (1-x_retain_list[i]) * (1-y_retain_list[i]) * A1_list[i] + \
                               x_retain_list[i]  * (1-y_retain_list[i]) * A2_list[i] + \
                            (1-x_retain_list[i]) *    y_retain_list[i]  * A3_list[i] + \
                               x_retain_list[i]  *    y_retain_list[i]  * A4_list[i] ;
        }

}




void generating_magnification_lightcurve(float *hjds,float *lightcurve,register int nhjd,\
                                        float *parmfit,\
                                        struct whether_offset_4mag *one_map_all_layer_struct,\
                                        /*float *all_layer_corner_mag,\
                                        _Bool *all_layer_whether_densed,\
                                        int *all_layer_sequence_number_in_next_layer_file,\*/
                                        int *layer_length_accumulated_in_front_of_each_layer,\
                                        float box_size, float s, float q)
{
        float t0,u0,te,alpha;
        t0 = parmfit[0];
        u0 = parmfit[1];
        te = parmfit[2];
        alpha = parmfit[3];

        float cos_alpha,sin_alpha;
        cos_alpha = cos(alpha/180*3.1415926);
        sin_alpha = sin(alpha/180*3.1415926);

        float sum_x_list[nhjd] ; 
        float sum_y_list[nhjd] ; 
        for(register int i=0;i<nhjd;i++) // vectorize
        {
                sum_x_list[i] = cos_alpha * (hjds[i]-t0)/te - u0*sin_alpha ;
                sum_y_list[i] = sin_alpha * (hjds[i]-t0)/te + u0*cos_alpha ;
        }


        float s_square ;
        _Bool whether_in_map[nhjd] ;
        for(register int i=0;i<nhjd;i++) // vectorize
        {
                // set the initial value to be true
                whether_in_map[i] = 1 ;
        }
        register int n_in_map = nhjd ;
        for(register int i=0;i<nhjd;i++)
        {
                // boxsize, sum_x_list[i] and sum_y_list[i] are read from cache
                if( sum_x_list[i]<=-box_size || sum_x_list[i]>=box_size || sum_y_list[i]<=-box_size || sum_y_list[i]>=box_size )
                {
                        //s = pow(sum_x_list[i]*sum_x_list[i]+sum_y_list[i]*sum_y_list[i],0.5) ;
                        //lightcurve[i] = (s*s+2) / s / pow(s*s + 4,0.5) ;
                        s_square = sum_x_list[i]*sum_x_list[i]+sum_y_list[i]*sum_y_list[i] ;
                        if(s>4.25)
                        {
                                s_square *= (1.0+q) ;
                        }
                        lightcurve[i] = (s_square+2) / sqrt( s_square * (s_square + 4) ) ;

                        n_in_map -= 1 ;
                        whether_in_map[i] = 0 ;
                }
                
        }

        if( n_in_map>0 )
        {
                float in_map_x_list[n_in_map] ;
                float in_map_y_list[n_in_map] ;
                float in_map_A_list[n_in_map] ;
                register int j=0 ;
                for(register int i=0;i<nhjd;i++)
                {
                        if(whether_in_map[i])
                        {
                                in_map_x_list[j] = sum_x_list[i] ;
                                in_map_y_list[j] = sum_y_list[i] ;
                                j += 1 ;
                        }
                }




                /******* change here ******/
                register int Npoint = 8 ;
                //printf("%d\n",Npoint) ;
                /*********** end **********/

                register int quotient ;
                register int remainder ;

                quotient  = n_in_map / Npoint ;
                remainder = n_in_map % Npoint ;

                float x_list[Npoint] ; // doesn't matter if x_list is longer
                float y_list[Npoint] ; 
                float A_list[Npoint] ; 

                if(quotient>0)
                {
                        for(register int i=0;i<quotient;i++)
                        {       
                                for(register int j=0;j<Npoint;j++)
                                {
                                        x_list[j] = in_map_x_list[i*Npoint+j] ;
                                        y_list[j] = in_map_y_list[i*Npoint+j] ;
                                }
                                
                                interpolating_to_get_8_magnification(x_list,y_list,A_list,\
                                                        one_map_all_layer_struct,\
                                                        /*all_layer_corner_mag,\
                                                        all_layer_whether_densed,\
                                                        all_layer_sequence_number_in_next_layer_file,\*/
                                                        layer_length_accumulated_in_front_of_each_layer,\
                                                        box_size) ;
                                /*interpolating_to_get_N_magnification(x_list,y_list,A_list,\
                                                        Npoint,\
                                                        all_layer_corner_mag,\
                                                        all_layer_whether_densed,\
                                                        all_layer_sequence_number_in_next_layer_file,\
                                                        layer_length_accumulated_in_front_of_each_layer,\
                                                        box_size) ;*/

                                for(register int j=0;j<Npoint;j++)
                                {
                                        in_map_A_list[i*Npoint+j] = A_list[j] ;
                                }
                        }
                }
                
                if(remainder>0)
                {
                        for(register int j=0;j<remainder;j++)
                        {
                                x_list[j] = in_map_x_list[quotient*Npoint+j] ;
                                y_list[j] = in_map_y_list[quotient*Npoint+j] ;
                        }
                        
                        interpolating_to_get_N_magnification(x_list,y_list,A_list,\
                                                remainder,\
                                                one_map_all_layer_struct,\
                                                /*all_layer_corner_mag,\
                                                all_layer_whether_densed,\
                                                all_layer_sequence_number_in_next_layer_file,\*/
                                                layer_length_accumulated_in_front_of_each_layer,\
                                                box_size) ;

                        for(register int j=0;j<remainder;j++)
                        {
                                in_map_A_list[quotient*Npoint+j] = A_list[j] ;
                        }  
                }

                //register int j=0 ;
                j = 0 ;
                for(register int i=0;i<nhjd;i++)
                {
                        if( whether_in_map[i]==1 )
                        {
                                lightcurve[i] = in_map_A_list[j] ;
                                j += 1 ;
                        }
                }

        }

}


/*float interpolating_to_get_magnification(float x,float y,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size)
{
	if( x<=-box_size || x>=box_size || y<=-box_size || y>=box_size )
        {
            float s = pow(x*x+y*y,0.5);
            return (s*s+2) / s / pow(s*s + 4,0.5);
        }
 	
	float A ;
	
	x = 0.5 * ( x + box_size ) / box_size ;
	y = 0.5 * ( y + box_size ) / box_size ;

	if ( all_layer_whether_densed[0] )
	{
		bool arg_x = int( 2 * x ) ;
		bool arg_y = int( 2 * y ) ;

		int sequence_number_in_next_layer_file = all_layer_sequence_number_in_next_layer_file[ layer_length_accumulated_in_front_of_each_layer[0] + 0 * 4 + 0 ] ;

		A = checking_whether_densed_to_next_layer(x,y,arg_x,arg_y,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,sequence_number_in_next_layer_file,0+1) ;
	}
	else
	{
		float A1 , A2 , A3 , A4 ;
                A1 = all_layer_corner_mag[ 0 ] ;
                A2 = all_layer_corner_mag[ 1 ] ;
                A3 = all_layer_corner_mag[ 2 ] ;
                A4 = all_layer_corner_mag[ 3 ] ;

                A = (1-x)*(1-y)*A1 + x*(1-y)*A2 + (1-x)*y*A3 + x*y*A4 ;

	}
	
	return A ;
} */

/*void generating_magnification_lightcurve(float *hjds,float *lightcurve,float *xtraj_list,float *ytraj_list,int nhjd,float *parmfit,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size)
{
        float t0,u0,te,alpha;
        t0 = parmfit[0];
        u0 = parmfit[1];
        te = parmfit[2];
        alpha = parmfit[3];

        float cos_alpha,sin_alpha;
        cos_alpha = cos(alpha/180*3.1415926);
        sin_alpha = sin(alpha/180*3.1415926);

	float xtraj , ytraj ;

	for(int i=0;i<nhjd;i++)
        {
                xtraj = cos_alpha * (hjds[i]-t0)/te - u0*sin_alpha ;
                ytraj = sin_alpha * (hjds[i]-t0)/te + u0*cos_alpha ;
		
		xtraj_list[i] = xtraj;
		ytraj_list[i] = ytraj;

		lightcurve[i] = interpolating_to_get_magnification(xtraj,ytraj,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);
	}

}*/

float getchi2_sub(float *iflux, float *iferr, float *lc, int nhjd)
{
        double d0=0,d1=0,b00=0,b01=0,b10=0,b11=0;
        double wght,det;
        double fs,fb;
        for(int i=0;i<nhjd;i++)
        {
                wght = iflux[i]/iferr[i]/iferr[i];
                d0 += wght*lc[i];
                d1 += wght;
                b00 += lc[i]*lc[i]/iferr[i]/iferr[i];
                b01 += lc[i]/iferr[i]/iferr[i];
                b11 += 1/iferr[i]/iferr[i];
        }
        b10 = b01;
        //std::cout<<b00<<" "<<b01<<" "<<b11<<std::endl;

        fs = b11*d0-b01*d1;
        fb = -b10*d0+b00*d1;
        det = b00*b11-b10*b01;
        if(det != 0)
        {
                fs /= det;
                fb /= det;
        }
        else
        {
                return -1;
                //fs = 0;
                //fb = d0/b11;
        }
        float res,chi2=0;
        for(int i=0;i<nhjd;i++)
        {
                //std::cout<<(iflux[i]-fb)/fs<<" "<<lc[i]<<std::endl;
                res = iflux[i]-lc[i]*fs-fb;
                //cout << res << endl;
                chi2 += res*res/iferr[i]/iferr[i];
        }
	return chi2;
}

float getchi2(float *hjds,float *iflux, float *iferr,int nhjd,float *parmfit,\
                struct whether_offset_4mag *one_map_all_layer_struct,\
                /*float *all_layer_corner_mag,\
                _Bool *all_layer_whether_densed,\
                int *all_layer_sequence_number_in_next_layer_file,\*/
                int *layer_length_accumulated_in_front_of_each_layer,\
                float box_size, float s, float q)
{
        //clock_t start_clock, end_clock ;

        //start_clock = clock() ;

        //float *lightcurve = new float[nhjd];
        float lightcurve[nhjd] ;
        float chi2;

        generating_magnification_lightcurve(hjds,lightcurve,nhjd,parmfit,\
                                                one_map_all_layer_struct,\
                                                /*all_layer_corner_mag,\
                                                all_layer_whether_densed,\
                                                all_layer_sequence_number_in_next_layer_file,\*/
                                                layer_length_accumulated_in_front_of_each_layer,\
                                                box_size, s, q);

        chi2 = getchi2_sub(iflux,iferr,lightcurve,nhjd);
        
	//delete [] lightcurve;

        //end_clock = clock() ;
        //double clock_taken = (double)( end_clock - start_clock ) / CLOCKS_PER_SEC ;
        //printf("Npoint need %e second",clock_taken) ;


	return chi2;
}

// to do : 
//test small granularity time in C code; 
//use valgrind to test; 
//statistics on C code time for different maps in a whole grid search

/*
int main()
{

}*/

/*int printlc(float *hjds,int nhjd, float *parmfit,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size)
{
        float *lc = new float[nhjd];
        float *xtraj_list = new float[nhjd];
        float *ytraj_list = new float[nhjd];

        generating_magnification_lightcurve(hjds,lc,xtraj_list,ytraj_list,nhjd,parmfit,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);
        //std::ofstream outFile("lc.txt");
        //for(int i=0;i<nhjd;i++)
        //{
        //        outFile<<hjds[i]<<" "<<lc[i]<<" "<<xtraj_list[i]<<" "<<ytraj_list[i]<<" "<<parmfit[0]<<" "<<parmfit[1]<<" "<<parmfit[2]<<" "<<parmfit[3]<<" "<<std::endl;
        //}
        //outFile.close();
        delete [] lc;
        delete [] xtraj_list;
        delete [] ytraj_list;

	return 1;

}*/

//extern "C"
//{

void * wrapgetchi2(float *hjds,float *iflux, float *iferr,int nhjd,float *parmfit,\
                        struct whether_offset_4mag *one_map_all_layer_struct,\
                        /*float *all_layer_corner_mag,\
                        _Bool *all_layer_whether_densed,\
                        int *all_layer_sequence_number_in_next_layer_file,\*/
                        int *layer_length_accumulated_in_front_of_each_layer,\
                        float box_size,\
                        float *chi2, float s, float q)
{
        *chi2 = getchi2(hjds,iflux,iferr,nhjd,parmfit,\
                        one_map_all_layer_struct,\
                        /*all_layer_corner_mag,\
                        all_layer_whether_densed,\
                        all_layer_sequence_number_in_next_layer_file,\*/
                        layer_length_accumulated_in_front_of_each_layer,\
                        box_size, s, q);
}

/*void * wrapprintlc(float *hjds,int nhjd,float *parmfit,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size)
{
        printlc(hjds,nhjd,parmfit,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);
}

void * wrapinterpolating_to_get_magnification(float x,float y,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size,float *A)
{
        *A = interpolating_to_get_magnification(x,y,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);
}*/

//}











