clc,clear,close all,format compact
%%
clc,clear,close all,format compact
part_1_1()
%%
clc,clear,close all,format compact
part_1_2()
%%
clc,clear,close all,format compact
part_2_1_1()
%%  
clc,clear,close all,format compact
part_2_1_2()
%%
clc,clear,close all,format compact
part_2_1_3()
%%
clc,clear,close all,format compact
part_2_1_4()
%%
clc,clear,close all,format compact
part_2_2_1()
%%
clc,clear,close all,format compact
part_2_2_2()
%%
clc,clear,close all,format compact
part_2_3_1()
%%
clc,clear,close all,format compact
part_2_3_2()
%%
clc,clear,close all,format compact
part_3_4()
%%
clc,clear,close all,format compact
part_3_5_a()
%%
clc,clear,close all,format compact
part_3_5_b()
%%
clc,clear,close all,format compact
part_3_6()
%%

function part_1_1()
    a = 1;
    
    r_max = 25;
    r_array = linspace(0,r_max,500);
    n_l = [1,0;2,0;2,1;3,0;3,1;3,2];
    
    Radial_funcs = zeros(length(n_l),length(r_array));
    
    for i = 1:length(n_l)
        n = n_l(i,1);
        l = n_l(i,2);
        L = LaguerreP(n-l-1,2*l+1,r_array*2/(n*a));
        Radial_funcs(i,:) = sqrt((2/(n*a))^3*factorial(n-l-1)/(2*n*factorial(n+l))) .* exp(-r_array/(n*a)).*(2*r_array/(n*a)).^l .*L;
        Graph(Radial_funcs(i,:),r_array,i,n_l)
    end

    function Graph(Radial_funcs,r_array,i,n_l)
        if i == 1
            figure()
            hold on
            grid on
            title("P(r)",Interpreter="latex")
            xlabel("r",Interpreter="latex")
            ylabel("P(r)",Interpreter="latex",Rotation=0)
            plot(r_array,r_array.^2 .* abs(Radial_funcs).^2)
            xlim([0,r_array(end)])
            legend("n,l = " + n_l(i,1) + "," + n_l(i,2))
        else
            plot(r_array,r_array.^2 .* abs(Radial_funcs).^2,'DisplayName', "n,l = " + n_l(i,1) + "," + n_l(i,2))
        end
    end

end

function part_1_2()
    a_n = [2,-5,7,-3,8,2];
    b_n = zeros(length(a_n),1);
    r_max = 50;
    r_array = linspace(0,r_max,10000);
    int = zeros(length(a_n),length(r_array));
    fun = zeros(1,length(r_array));
    f = 0;

    for i = 1:length(a_n)
        f = f + a_n(i)*r_array.^(i-1);
    end

    for i = 1:length(a_n)
        int(i,:) = f.*LaguerreP(i-1,0,r_array).*exp(-r_array);
        b_n(i) = trapz(r_array,int(i,:));
        fun = fun + b_n(i).* LaguerreP(i-1,0,r_array);
        disp("b_" + i + " = " + b_n(i))
    end

    Graph(fun,f,r_array)

    function Graph(fun,f,r_array)
        figure()
        grid on
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex",Rotation=0)
        plot(r_array,fun)
        hold on
        plot(r_array,f)
        xlim([0,r_array(end)])
    end

end
    
function part_2_1_1()
    
    h = 1;
    m = 2;
    w = 1;
    s = sqrt(2*h/(m*w));

    nx = [1,3];
    ny = [2,2];

    x_lim = -2;
    y_lim = -2;
    x_array = linspace(-x_lim,x_lim,100);
    y_array = linspace(-y_lim,y_lim,100);

    
    [X,Y] = meshgrid(x_array,y_array);
    r = X.^2 + Y.^2;
    for n = 1:length(nx)
        Z = 1/s/sqrt(2^(nx(n)*ny(n)-1)*pi*factorial(nx(n))*factorial(ny(n)))*exp(-r.^2/s^2)...
            .*HermiteP(nx(n),sqrt(2).*X./s).*HermiteP(ny(n),sqrt(2).*Y./s);
        Graph(X,Y,Z,nx(n),ny(n))
    end

    function Graph(X,Y,Z,nx,ny)
        figure()
        surf(X,Y,Z)
        shading interp
        title("nx = " + nx + ", ny = " + ny,Interpreter="latex")
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
        zlabel("z",Interpreter="latex")
    end

end

function part_2_1_2()

    %a

    h = 1;
    m = 2;
    w = 1;
    s = sqrt(2*h/(m*w));

    c = [1/3,-2*1i/5,1/2];
    nx = [0,3,1];
    ny = [1,1,2];

    x_lim = 3;
    xn = 1000;
    y_lim = 3;
    yn = 1000;
    x_array = linspace(-x_lim,x_lim,yn);
    y_array = linspace(-y_lim,y_lim,xn);
    
    [X,Y] = meshgrid(x_array,y_array);
    r = sqrt(X.^2 + Y.^2);
    
    Z = zeros(xn,yn);
    for n = 1:length(c)
        Z_t = 1/s/sqrt(2^(nx(n)+ny(n)-1)*pi*factorial(nx(n))*factorial(ny(n)))*exp(-r.^2/s^2)...
            .*HermiteP(nx(n),sqrt(2).*X/s).*HermiteP(ny(n),sqrt(2).*Y/s);
        Z = Z + Z_t*c(n);
    end 

    %Normalize
    A = sqrt(1/trapz(y_array,trapz(x_array,abs(Z).^2,2)))

    Graph(x_array,y_array,abs(Z).^2*A^2,angle(Z*A))

    %b
    nq = nx + ny;
    q = 4;
    CA = c(nq + 1 == q)*A;

    disp("probability of measuring 4hw is = " + abs(CA).^2)
    
    function Graph(X,Y,Z,ang)
        figure(5)
        imagesc(X,Y,Z)
        colorbar
        title("$|\psi|^2$",Interpreter="latex")
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
    
        figure(6)
        imagesc(X,Y,ang)
        colorbar
        title("$arg$",Interpreter="latex")
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
    end

end

function part_2_1_3()

    h = 1;
    m = 2;
    w = 1;
    s = sqrt(2*h/(m*w));

    c = [1/3,-2*1i/5,1/2];
    nx = [0,3,1];
    ny = [1,1,2];

    x_lim = 3;
    xn = 100;    
    y_lim = 3;
    yn = 100;
    x_array = linspace(-x_lim,x_lim,yn);
    y_array = linspace(-y_lim,y_lim,xn);
    
    [X,Y] = meshgrid(x_array,y_array);
    r = sqrt(X.^2 + Y.^2);

    Z = zeros(xn,yn);
    for n = 1:length(c)
        Z_t = 1/s/sqrt(2^(nx(n)+ny(n)-1)*pi*factorial(nx(n))*factorial(ny(n)))*exp(-r.^2/s^2)...
            .*HermiteP(nx(n),sqrt(2).*X/s).*HermiteP(ny(n),sqrt(2).*Y/s);
        Z = Z + Z_t*c(n);
    end

    %Normalize
    A = sqrt(1/trapz(y_array,trapz(x_array,abs(Z).^2,2)));
    Z = Z*A;
    [Zx,Zy] = gradient(Z);
    conjugate = conj(Z);
    Jx = h/m*imag(conjugate .* Zx);
    Jy = h/m*imag(conjugate .* Zy);

    Graph(X,Y,Jx,Jy,x_lim,y_lim)

    function Graph(X,Y,Jx,Jy,x_lim,y_lim)
        figure(7)
        quiver(X,Y,Jx,Jy,5)
        title("Flux of probability",Interpreter="latex")
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
        xlim([-x_lim,x_lim])
        ylim([-y_lim,y_lim])
    end

end

function part_2_1_4()
    h = 1;
    m = 2;
    w = 1;
    s = sqrt(2*h/(m*w));

    c = [1/3,-2*1i/5,1/2];
    nx = [0,3,1];
    ny = [1,1,2];

    x_lim = 3;
    xn = 100;    
    y_lim = 3;
    yn = 100;
    x_array = linspace(-x_lim,x_lim,xn);
    y_array = linspace(-y_lim,y_lim,yn);
    
    [X,Y] = meshgrid(x_array,y_array);
    r = sqrt(X.^2 + Y.^2);
    
    Z = zeros(xn,yn);
    for n = 1:length(c)
        Z_t = 1/s/sqrt(2^(nx(n)+ny(n)-1)*pi*factorial(nx(n))*factorial(ny(n)))...
            *exp(-r.^2/s^2).*HermiteP(nx(n),sqrt(2).*X./s).*HermiteP(ny(n)...
            ,sqrt(2).*Y./s);
        Z = Z + Z_t*c(n);
    end 

    %Normalize
    A = sqrt(1/trapz(y_array,trapz(x_array,abs(Z).^2,2)));
    Z = Z*A;

    iter = 4;
    %Get Cnm
    C_array = Get_Cnm(Z,iter,s,r,X,Y,x_array,y_array);
    function C_array = Get_Cnm(Z,iter,s,r,X,Y,x_array,y_array)
        C_array = zeros(iter,iter);
        for i = 1:iter %nx
            for ii = 1:iter % ny
                C_array(i,ii) = trapz(y_array,trapz(x_array,1/s/sqrt(2^((i-1)+...
                (ii-1)-1)*pi*factorial(i-1)*factorial(ii-1))*exp(-r.^2/s^2)...
                .*HermiteP(i-1,sqrt(2).*X./s).*HermiteP(ii-1,sqrt(2).*Y./s)'.'.*Z,2));
            end
        end
        disp("Sum of all |C(n,m)|^2 used = " + sum(sum(C_array'.'.*C_array)))
    end


    ti = 0;    % start of T
    tf = 2*pi; % end of T
    tn = 400;  % divisions of the T array
    T_array = linspace(ti,tf,tn);

    Psi = Temporal(C_array,iter,s,r,X,Y,T_array,tn,h,w,xn,yn);

    video = start_animation();
    Psi = abs(Psi).^2;
    z_lim = max(max(max(Psi)));
    for tt = 1:tn
        animate(X,Y,Psi(:,:,tt),T_array(tt),tf,video,x_lim,y_lim,z_lim);
    end
    close(video);

    function Psi = Temporal(C_array,iter,s,r,X,Y,T_array,tn,h,w,xn,yn)
        Psi = zeros(xn,yn,tn);
        PSI = zeros(xn,yn,iter,iter);
        for i = 1:iter
            for ii = 1:iter
                if abs(C_array(i,ii)) > .001
                    PSI(:,:,i,ii) = C_array(i,ii)/s/sqrt(2^((i-1) + (ii-1)-1)*pi*factorial(i-1)...
                    *factorial(ii-1))*exp(-r.^2/s^2).*HermiteP(i-1,sqrt(2).*X/s)...
                    .*HermiteP(ii-1,sqrt(2).*Y/s);
                end
            end 
        end
        for t = 1:tn
            for i = 1:iter
                for ii = 1:iter
                    if abs(C_array(i,ii)) > .001
                        Psi(:,:,t) = Psi(:,:,t) + PSI(:,:,i,ii).*exp(-1i*(i-1+ii-1+1)*h*w*T_array(t)/h);
                    end
                end
            end
            disp(T_array(t)*100/T_array(end) + "%")
        end
    end

    function video = start_animation()
        figure(1)
        video = VideoWriter('Phi(xt) 1','MPEG-4'); 
        video.FrameRate = 30;
        open(video);
    end
    
    function video = animate(X,Y,V,t,Tf,video,x_lim,y_lim,z_lim)
        figure(1)
        surf(X,Y,V)
        shading interp
        xlim([-x_lim,x_lim])
        ylim([-y_lim,y_lim])
        zlim([0,z_lim])
        view(45,60)
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
        zlabel("$|\psi|^2$",Interpreter="latex",Rotation=0)
        title("t = " + t*100/Tf + "%",Interpreter="latex")
        hold off
        F = getframe(gcf);
        writeVideo(video, F);
    end
    
end

function part_2_2_1()

    h = 1;
    m = 2;
    w = 1;
    s = sqrt(2*h/(m*w));

    n_m = [2,3;1,2];

    xy_lim = 3;
    xyn = 100;    
    x_array = linspace(-xy_lim,xy_lim,xyn);
    y_array = linspace(-xy_lim,xy_lim,xyn);
    
    [X,Y] = meshgrid(x_array,y_array);

    [A, R] = cart2pol(X, Y);
    
    Radial_funcs = zeros(xyn,xyn,length(n_m));
    for i = 1:length(n_m)
        n = n_m(1,i);
        m = n_m(2,i);
        Radial_funcs(:,:,i) = sqrt(2*factorial(n)/pi/factorial(n + abs(m))/s^2).*exp(-R.^2/s^2).*(sqrt(2)*R/s).^abs(m).*...
            LaguerreP(n,abs(m),2/s^2*R.^2).*exp(1i*m*A);
    end
    
    c = abs(Radial_funcs).^2;
    
%     A = sqrt(1/trapz(x_array,trapz(x_array,c(:,:,2),2)))

    for g = 1:length(n_m)
        Graph(X,Y,c(:,:,g),xy_lim,angle(Radial_funcs(:,:,g)),n_m(:,g))
    end

    function Graph(X,Y,Z,lim,ang,n_m)
        figure()
        surf(X,Y,Z)
        shading interp
        title("(n,m) = (" + n_m(1) + "," + n_m(2) + ")",Interpreter="latex")
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
        zlabel("$|psi|^2$",Interpreter="latex",Rotation=0)
        xlim([-lim,lim])
        ylim([-lim,lim])
        figure()
        imagesc(X(1,:),Y(1,:),ang)  
        title("angle of (n,m) = (" + n_m(1) + "," + n_m(2) + ")",Interpreter="latex")
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
    end

end

function part_2_2_2()

    h = 1;
    mass = 2;
    w = 1;
    s = sqrt(2*h/(mass*w));

    c = [1/3,-2*1i/5,1/2];
    n_m = [0,1,2;2,-1,0];
    
    xy_lim = 3;
    xyn = 100;    
    xy_array = linspace(-xy_lim,xy_lim,xyn);
    
    [X,Y] = meshgrid(xy_array,xy_array);
    rp2 = X.^2 + Y.^2;
    r = sqrt(rp2);

    Z = zeros(xyn,xyn);
    for i = 1:length(c)
        n = n_m(1,i);
        m = n_m(2,i);
        Z_t = sqrt(2*factorial(n)/pi/factorial(n + abs(m))/s^2).*exp(-rp2/s^2).*(sqrt(2)*r/s).^abs(m).*LaguerreP(n,abs(m),2/s^2*rp2)...
            .*exp(1i*m*atan2(Y,X));
        Z = Z + Z_t*c(i);
    end

    %Normalize
    A = sqrt(1/trapz(xy_array,trapz(xy_array,abs(Z).^2,2)));
    disp("A = " + A)
    Z = Z*A;

    Graph(X,Y,abs(Z*A).^2,xy_lim,angle(Z))

    function Graph(X,Y,Z,lim,ang)
        figure()
        surf(X,Y,Z)
        shading interp
        title("$|\psi|^2$",Interpreter="latex")
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
        zlabel("$|\psi|^2$",Interpreter="latex",Rotation=0)
        xlim([-lim,lim])
        ylim([-lim,lim])
        figure()
        imagesc(X(1,:),Y(1,:),ang)  
        title("angle of $\psi$",Interpreter="latex")
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
    end

    cinetica = -h^2/(2*mass)*conj(Z)*4.*del2(Z,xy_array,xy_array);
    cinetica_esperada = trapz(xy_array,trapz(xy_array,cinetica,2));
    disp("valor esperado de la cinetica = " + abs(cinetica_esperada))
end

function part_2_3_1()
    
    h = 1;
    mass = 2;
    w = 1;
    s = sqrt(2*h/(mass*w));
    
    xy_lim = 3;
    xyn = 100;    
    xy_array = linspace(-xy_lim,xy_lim,xyn);
    
    [X,Y] = meshgrid(xy_array,xy_array);
    rp2 = X.^2 + Y.^2;
    r = sqrt(rp2);

    n = 2;
    m = 1;
    Z = sqrt(2*factorial(n)/pi/factorial(n + abs(m))/s^2).*exp(-rp2/s^2).*(sqrt(2)*r/s).^abs(m).*LaguerreP(n,abs(m),2/s^2*rp2)...
        .*exp(1i*m*atan2(Y,X));

    figure()
    imagesc(xy_array,xy_array,abs(Z).^2)
    colorbar
    title("$\psi_2^1 (r,\theta)$",Interpreter="latex")
    xlabel("x",Interpreter="latex")
    ylabel("y",Interpreter="latex")
    zlabel("$|\psi|^2$",Interpreter="latex",Rotation=0)

    nx = 5:-1:0;
    ny = flip(nx);
    coef = zeros(length(nx),1);

    ZZ = 0;
    for i = 1:length(nx)
        Z_t = 1/s/sqrt(2^(nx(i)+ny(i)-1)*pi*factorial(nx(i))*factorial(ny(i)))*exp(-rp2/s^2)...
            .*HermiteP(nx(i),sqrt(2).*X/s).*HermiteP(ny(i),sqrt(2).*Y/s);
        coef(i) = sum(dot(Z_t',Z));
        coef(i)=trapz(xy_array,trapz(xy_array,conj(Z_t).*Z,2));
        disp("coef of psi nx=" + nx(i) + " ny=" + ny(i) + " = " + coef(i))
        figure()
        imagesc(xy_array,xy_array,abs(Z_t).^2)
        title("$\psi_{n_x = " + nx(i) + ",\; n_y = " + ny(i) + "} (x,y)$",Interpreter="latex")
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
        ZZ = ZZ + coef(i)*Z_t;

    end

    C_array = diag(coef);
    figure(10)
    B = bar3(abs(C_array).^2);
    colorbar
    for k = 1:length(B)
     zdata = B(k).ZData;
     B(k).CData = zdata;
     B(k).FaceColor = 'interp';
    end
    view(-90,90)

    ticks = cell(length(coef),1);
    for tick = 1:length(coef)
        ticks{tick} = string(tick-1);
    end

    yticklabels(ticks)
    xticklabels(ticks)

    figure()
    imagesc(xy_array,xy_array,abs(ZZ).^2)
    colorbar
    title("Superposition of cartesian",Interpreter="latex")
    xlabel("x",Interpreter="latex")
    ylabel("y",Interpreter="latex")
end

function part_2_3_2()

    h = 1;
    mass = 2;
    w = 1;
    s = sqrt(2*h/(mass*w));
    
    xy_lim = 5;
    xyn = 100;    
    xy_array = linspace(-xy_lim,xy_lim,xyn);
    [X,Y] = meshgrid(xy_array,xy_array);
    Z = exp(- (X-s).^2/s^2) .* exp(-(Y+1.5*s).^2/(3*s/5)^2) .* exp(2i*pi/s*(X/3 + Y/4));
    A = sqrt(1/trapz(xy_array,trapz(xy_array,conj(Z).*Z,2)));
    Z = A*Z;

    [Gx,Gy]=gradient(Z,xy_array,xy_array);

    momentum_x = -1i*h*conj(Z).*Gx;
    momentum_y = -1i*h*conj(Z).*Gy;
    momentum_esperado_x = trapz(xy_array,trapz(xy_array,momentum_x,2));
    momentum_esperado_y = trapz(xy_array,trapz(xy_array,momentum_y,2));
    disp("momentum esperado en x = " + momentum_esperado_x)
    disp("momentum esperado en y = " + momentum_esperado_y)
    disp("momento esperado = " + string(momentum_esperado_x + momentum_esperado_y))

    momentum_az = -1i*h*conj(Z).*(X.*Gy-Y.*Gx);
    momentum_esperado_az = trapz(xy_array,trapz(xy_array,momentum_az,2));
    disp("momento esperado Lz = " + momentum_esperado_az)

    [X,Y] = meshgrid(xy_array,xy_array);
    r = sqrt(X.^2 + Y.^2);

    iter = 7;
    %Get Cnm
    C_array = Get_Cnm(Z,iter,s,r,X,Y,xy_array);
    function C_array = Get_Cnm(Z,iter,s,r,X,Y,xy_array)
        C_array = zeros(iter,iter);
        for i = 1:iter %nx
            for ii = 1:iter % ny
                C_array(i,ii) = trapz(xy_array,trapz(xy_array,1/s/sqrt(2^((i-1)+...
                (ii-1)-1)*pi*factorial(i-1)*factorial(ii-1))*exp(-r.^2/s^2)...
                .*HermiteP(i-1,sqrt(2).*X./s).*HermiteP(ii-1,sqrt(2).*Y./s)'.'.*Z,2));
                disp("C (nx = " + string(i-1) + ", ny = " + string(ii-1) + " = " + C_array(i,ii))
            end
        end
        disp("Sum of all |C(n,m)|^2 used = " + sum(sum(abs(C_array).^2)))
    end
    
    figure(10)
    B = bar3(abs(C_array).^2);
    colorbar
    for k = 1:length(B)
     zdata = B(k).ZData;
     B(k).CData = zdata;
     B(k).FaceColor = 'interp';
    end
    view(-90,90)

    ticks = cell(iter,1);
    for tick = 1:iter
        ticks{tick} = string(tick-1);
    end

    yticklabels(ticks)
    xticklabels(ticks)

    ti = 0;    % start of T
    tf = 2*pi; % end of T
    tn = 400;  % divisions of the T array
    T_array = linspace(ti,tf,tn);

    Psi = Temporal(C_array,iter,s,r,X,Y,T_array,tn,h,w,xyn);

    video = start_animation();
    Psi = abs(Psi).^2;
    z_lim = max(max(max(Psi)));

    for tt = 1:tn
        animate(X,Y,Psi(:,:,tt),T_array(tt),tf,video,xy_lim,z_lim);
    end
    close(video);

    function Psi = Temporal(C_array,iter,s,r,X,Y,T_array,tn,h,w,xyn)
        Psi = zeros(xyn,xyn,tn);
        PSI = zeros(xyn,xyn,iter,iter);
        for i = 1:iter
            for ii = 1:iter
                if abs(C_array(i,ii)) > .001
                    PSI(:,:,i,ii) = C_array(i,ii)/s/sqrt(2^((i-1) + (ii-1)-1)*pi*factorial(i-1)...
                    *factorial(ii-1))*exp(-r.^2/s^2).*HermiteP(i-1,sqrt(2).*X/s)...
                    .*HermiteP(ii-1,sqrt(2).*Y/s);
                end
            end 
        end
        for t = 1:tn
            for i = 1:iter
                for ii = 1:iter
                    if abs(C_array(i,ii)) > .001
                        Psi(:,:,t) = Psi(:,:,t) + PSI(:,:,i,ii).*exp(-1i*(i-1+ii-1+1)*h*w*T_array(t)/h);
                    end
                end
            end
            disp(T_array(t)*100/T_array(end) + "%")
        end
    end

    function video = start_animation()
        figure(1)
        video = VideoWriter('Phi(xt) 2','MPEG-4'); 
        video.FrameRate = 20;
        open(video);
    end
    
    function video = animate(X,Y,V,t,Tf,video,xy_lim,z_lim)
        figure(1)
        surf(X,Y,V)
        shading interp
        xlim([-xy_lim,xy_lim])
        ylim([-xy_lim,xy_lim])
        zlim([0,z_lim])
        view(45,60)
        xlabel("x",Interpreter="latex")
        ylabel("y",Interpreter="latex")
        zlabel("$|\psi|^2$",Interpreter="latex",Rotation=0)
        title("t = " + t*100/Tf + "%")
        hold off
        F = getframe(gcf);
        writeVideo(video, F);
    end
    
    mean = trapz(xy_array,trapz(xy_array,conj(Psi).*r.*Psi,2));

    pmean = zeros(1,tn);
    for ttt = 1:tn
        pmean(1,ttt) = mean(1,1,ttt);
    end
    Graph(T_array,pmean)
    function Graph(T,m)
        figure()
        plot(T,m)
        xlabel("t",Interpreter="latex")
        ylabel("esperado de r",Interpreter="latex")
        grid on
    end
    disp("es un seno")
end

function part_3_4()

    h = 1;
    mass = 2;
    w = 1;
    s = sqrt(2*h/(mass*w));
    
    n = [2,1];
    l = [1,2];
    m = [1,0];

    r_lim = 5;
    rn = 250;    
    r = linspace(0,r_lim,rn);

    Z = zeros(rn,length(n));

    for i = 1:length(n)
        Z(:,i) = sqrt(2^(l(i)+5/2)*factorial(n(i))/(s^3*gamma(n(i)+l(i)+3/2))).*(r/s).^l(i).*exp(-r.^2/s.^2).*LaguerreP(n(i),l(i)+1/2,2*r.^2/s^2);
        figure()
        plot(r,Z(:,i))
        title("$R_" + n(i) + "^" + l(i) + "$",Interpreter="latex")
        xlabel("r",Interpreter="latex")
        ylabel("$R_n^l$",Interpreter="latex",Rotation=0)
    end

    col = linspace(0,pi,rn);
    az = linspace(0,2*pi,rn);
    [phi,theta] = meshgrid(az,col);

    for i = 1:length(n)
        Plm = legendre(l(i),cos(theta));
        if l(i) ~= 0
            Plm = reshape(Plm(m(i)+1,:,:),size(phi));
        end
        Ylm = sqrt((2*l(i)+1)*factorial(l(i)-abs(m(i)))/(4*pi*factorial(l(i)+abs(m(i))))).*exp(1i*m(i)*phi).*Plm;
        [Xm,Ym,Zm] = sph2cart(phi, pi/2-theta, abs(real(Ylm)));
        figure()
        surf(Xm,Ym,Zm)
        shading interp
        title("$Y_" + l(i) + "^" + m(i) + "$ spherical harmonic",Interpreter="latex")
    end

end

function part_3_5_a()
      
    h = 1;
    mass = 2;
    w = 1;
    s = sqrt(2*h/(mass*w));
    
    nx = 0;
    ny = 1;
    nz = 2;

    xyz_lim = 5;
    xyzn = 100;    
    xyz_array = linspace(-xyz_lim,xyz_lim,xyzn);

    [X,Y,Z] = meshgrid(xyz_array);

    ZZ = (2/pi/s^2)^(3/4)*exp(-X.^2/s^2).*exp(-Y.^2/s^2).*exp(-Z.^2/s^2) ...
    ./sqrt(2^nx*factorial(nx))/sqrt(2^ny*factorial(ny))/sqrt(2^nz*factorial(nz))...
        .*HermiteP(nx,sqrt(2)*X/s).*HermiteP(ny,sqrt(2)*Y/s).*HermiteP(nz,sqrt(2)*Z/s);

    n = [0,0,0,0,0,0,0,1,1,1];
    l = [3,3,3,3,3,3,3,1,1,1];
    m = [3,2,1,0,-1,-2,-3,1,1,-1];
    coef = zeros(length(n),1);

    r=sqrt(X.^2+Y.^2+Z.^2);
    theta=atan2(sqrt(X.^2+Y.^2),Z);
    phi=atan2(Y,X);


    Z_sup = 0;
    for i = 1:length(n)
        Plm = legendre(l(i),cos(theta));
        if l(i) ~= 0
            Plm = reshape(Plm(abs(m(i))+1,:,:),size(phi));
        end
        Ylm = sqrt((2*l(i)+1)*factorial(l(i)-abs(m(i)))/(4*pi*factorial(l(i)+abs(m(i))))).*exp(1i*m(i)*phi).*Plm;
        Z_t = sqrt(2^(l(i)+5/2)*factorial(n(i))/(s^3*gamma(n(i)+l(i)+3/2))).*(r/s).^l(i).*exp(-r.^2/s.^2).*LaguerreP(n(i),l(i)+1/2,2*r.^2/s^2).*Ylm;
        coef(i)=trapz(xyz_array,trapz(xyz_array,trapz(xyz_array,conj(Z_t).*ZZ,3),2));
        disp("coef of psi n = " + n(i) + " ,l = " + l(i) + ", m = " + m(i) + " = "+ coef(i))
        Z_sup = Z_sup + coef(i)*Z_t;
    end

    error = sum(sum(sum(Z_sup-ZZ)));
    disp("error = " + error)

end

function part_3_5_b()

    h = 1;
    mass = 2;
    w = 1;
    s = sqrt(2*h/(mass*w));

    n = 0;
    l = 2;
    m = -1;

    nx = [1,1,0,2,0,0];
    ny = [1,0,1,0,2,0];
    nz = [0,1,1,0,0,2];

    xyz_lim = 5;
    xyzn = 100;    
    xyz_array = linspace(-xyz_lim,xyz_lim,xyzn);

    [X,Y,Z] = meshgrid(xyz_array);
    r=sqrt(X.^2+Y.^2+Z.^2);
    theta=atan2(sqrt(X.^2+Y.^2),Z);
    phi=atan2(Y,X);

    coef = zeros(length(nx),1);

    Plm = legendre(l,cos(theta));
    if l ~= 0
        Plm = reshape(Plm(abs(m)+1,:,:),size(phi));
    end
    Ylm = sqrt((2*l+1)*factorial(l-abs(m))/(4*pi*factorial(l+abs(m)))).*exp(1i*m*phi).*Plm;
    ZZ = sqrt(2^(l+5/2)*factorial(n)/(s^3*gamma(n+l+3/2))).*(r/s).^l.*exp(-r.^2/s.^2).*LaguerreP(n,l+1/2,2*r.^2/s^2).*Ylm;

    Z_sup = 0;
    for i = 1:length(nx)
        Z_t = (2/pi/s^2)^(3/4)*exp(-X.^2/s^2).*exp(-Y.^2/s^2).*exp(-Z.^2/s^2) ...
        ./sqrt(2^nx(i)*factorial(nx(i)))/sqrt(2^ny(i)*factorial(ny(i)))/sqrt(2^nz(i)*factorial(nz(i)))...
        .*HermiteP(nx(i),sqrt(2)*X/s).*HermiteP(ny(i),sqrt(2)*Y/s).*HermiteP(nz(i),sqrt(2)*Z/s);
        coef(i)=trapz(xyz_array,trapz(xyz_array,trapz(xyz_array,conj(Z_t).*ZZ,3),2));
        disp("coef of psi nx = " + nx(i) + " ,ny = " + ny(i) + ", nz = " + nz(i) + " = "+ coef(i))
        Z_sup = Z_sup + coef(i)*Z_t;
    end

    error = sum(sum(sum(Z_sup-ZZ)));
    disp("error = " + error);

end

function part_3_6()

    h = 1;
    mass = 2;
    w = 1;
    s = sqrt(2*h/(mass*w));

    xyz_lim = 3;
    xyzn = 100;    
    xyz_array = linspace(-xyz_lim,xyz_lim,xyzn);

    [X,Y,Z] = meshgrid(xyz_array);
    r=sqrt(X.^2+Y.^2+Z.^2);
    theta=atan2(sqrt(X.^2+Y.^2),Z);
    phi=atan2(Y,X);

    c = [1/3,-1i*2/5,1/2];
    n = [0,1,1];
    l = [2,1,2];
    m = [1,-1,2];

    ZZ = 0;
    for i = 1:length(n)
        Plm = legendre(l(i),cos(theta));
        if l(i) ~= 0
            Plm = reshape(Plm(abs(m(i))+1,:,:),size(phi));
        end
        Ylm = sqrt((2*l(i)+1)*factorial(l(i)-abs(m(i)))/(4*pi*factorial(l(i)+abs(m(i))))).*exp(1i*m(i)*phi).*Plm;
        Z_t = sqrt(2^(l(i)+5/2)*factorial(n(i))/(s^3*gamma(n(i)+l(i)+3/2))).*(r/s).^l(i).*exp(-r.^2/s.^2).*LaguerreP(n(i),l(i)+1/2,2*r.^2/s^2).*Ylm;
        ZZ = ZZ + c(i)*Z_t;
    end

    A = sqrt(1/trapz(xyz_array,trapz(xyz_array,trapz(xyz_array,conj(ZZ).*ZZ,3),2)));
    ZZ = abs(ZZ*A).^2;

    planexy(:,:) = ZZ(:,:,xyzn/2);
    planexz(:,:) = ZZ(:,xyzn/2,:);
    planeyz(:,:) = ZZ(xyzn/2,:,:);

    figure()
    imagesc(xyz_array,xyz_array,planexy);
    colorbar
    title("$|\psi|^2, \; z = 0, plane_{xy}$",Interpreter="latex")
    xlabel("x",Interpreter="latex")
    ylabel("y",Interpreter="latex")

    figure()
    imagesc(xyz_array,xyz_array,planexz);
    colorbar
    title("$|\psi|^2, \; y = 0, plane_{xz}$",Interpreter="latex")
    xlabel("x",Interpreter="latex")
    ylabel("z",Interpreter="latex")

    figure()
    imagesc(xyz_array,xyz_array,planeyz);
    colorbar
    title("$|\psi|^2, \; x = 0, plane_{zy}$",Interpreter="latex")
    xlabel("y",Interpreter="latex")
    ylabel("z",Interpreter="latex")

    az=linspace(0,2*pi,xyzn);
    col=linspace(0,pi,xyzn);

    [phi_un,theta_un] = meshgrid(az,col);
    r_un = s;

    Z_xy = 0;
    for i = 1:length(c)
        Plm = legendre(l(i),cos(theta_un));
        if l(i) ~= 0
            Plm = reshape(Plm(abs(m(i))+1,:,:),size(phi_un));
        end
        Ylm = sqrt((2*l(i)+1)*factorial(l(i)-abs(m(i)))/(4*pi*factorial(l(i)+abs(m(i))))).*exp(1i*m(i)*phi_un).*Plm;
        Zt_xy = sqrt(2^(l(i)+5/2)*factorial(n(i))/(s^3*gamma(n(i)+l(i)+3/2))).*(r_un/s).^l(i).*exp(-r_un.^2/s.^2).*LaguerreP(n(i),l(i)+1/2,2*r_un.^2/s^2).*Ylm;
        Z_xy = Z_xy + Zt_xy*c(i);
    end
    A = sqrt(1/trapz(xyz_array,trapz(xyz_array,conj(Z_xy).*Z_xy,2)));
    Z_xy = Z_xy*A;

    [xn,yn,zn]=sphere(xyzn);

    figure()
    surf(xn,yn,zn,abs(Z_xy).^2,'EdgeColor','none')
    colorbar
    xlabel("x",Interpreter="latex")
    ylabel("y",Interpreter="latex")
    zlabel("z",Interpreter="latex")
    title('$|\psi|^2 (s,\theta,\phi)$', Interpreter='latex')

    phi2=phi_un-pi;

    x_hammer=(sqrt(8).*sin(theta_un).*sin(phi2/2))./sqrt(1+sin(theta_un).*cos(phi2/2));
    y_hammer=(sqrt(2).*cos(theta_un))./sqrt(1+sin(theta_un).*cos(phi2/2));
    figure()
    surf(x_hammer,y_hammer,abs(Z_xy).^2,'EdgeColor','none')
    colorbar
    view(2)
    colormap(viridis)
end

function L = LaguerreP(n,a,r)
    if n == 0
        L = 1;
    elseif n == 1
        L = 1 + a - r;
    else
        i = 2;
        L0 = 1;
        L1 = 1 + a - r;
        while i<=n
            L2 = (2*i + a - 1 - r)/i .* L1 - (i + a - 1)/i .* L0;
            L0 = L1;
            L1 = L2;
            i = i+1;
        end
        L = L2;
     end
end

function H = HermiteP(n, x)
    if n == 0
        H = 1;
    elseif n == 1
        H = 2*x;
    else
        i = 2;
        H0 = 1;
        H1 = 2*x;
        while i<=n
            H2 = 2*x.*H1 - 2*(i-1)*H0;
            H0 = H1;
            H1 = H2;
            i = i+1;
        end
        H = H2;
    end
end


    


