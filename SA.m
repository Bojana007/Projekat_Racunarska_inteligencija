
clear all; clc; close all;

%zadato nam je
T = 350;
Tc = 0.993;
Tmin = 0.1;

n = 2;
x_d = 0;
x_g = 5;

for i = 1:n
    curr(i) = x_d + ( x_g - x_d ) .* rand();
end

[result] = fitness_function( curr, 2 );

vec = [];

for i = 1:10

    T = 350;
    Tc = 0.993;
    Tmin = 0.1;

    [ curr, output, vektor ] = SA_4( curr, T, Tc, Tmin, n )

    if output < -1.85
        greska = curr
    end

    vec = [vec; output];
end

figure(1)
plot(vec)
yline(-1.8, 'r--') %minimum fje je -1.8
hold on;
figure(2)
boxplot( vec )


minimalno = min( vec )
maximalno = max( vec )
srednja_vr = mean( vec )
sd = std( vec )


function[ result ] = fitness_function( curr, n )
    
    result = 0;

    for i = 1:n
        result = result + (( sin( curr(i) ) * ( sin( ( i*curr(i)^2 ) / pi ) )^20 ) );
    end

    result = -result;
end


function[ curr, output, vektor ] = SA_4( curr, T, Tc, Tmin, n )
     
    vektor = [];

    while( T > Tmin)
        
        for i = 1:n
            if curr(i) < 0  
                curr(i) = 0;
            elseif curr(i) > 5 
                curr(i) = 5
            end
                
        end

       
        %treba da generisemo novi par
        for i = 1:n
            temp(i) = curr(i) + randn;
        end    

        %fitnes za privremeno resenje
        [ f_temp ] = fitness_function( temp, n );
        [ f_curr ] = fitness_function( curr, n );
       
        d_E =  f_temp - f_curr;
    
        if d_E < 0 
            
            curr = temp;
    
        else
    
            r = rand;
    
            if r < exp( - d_E / T )
                
                curr = temp;
                exp( - d_E / T )
                 
            end
    
        end
    
        T = Tc * T;
    
        [ output ] = fitness_function( curr, n );

       % vektor = [vektor, output];

        
    end
    
end







