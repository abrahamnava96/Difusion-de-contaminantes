function linspace(a,b,n)
    range(a,stop=b,length=n)
end

#rodear(A,a), rodea la matriz A de números a
function rodear(A,a)
    M,N = size(A)
    B = a*ones(M+2,N+2)
    B[2:M+1,2:N+1] = A
    return B
end

#id_normal(), identifica el vector normal exterior de un punto Ω_ij que esta en ∂Ω,
#de acuerdo a la clasificación de nodos de Ω

function id_normal(Ω,i,j)
    a = [Ω[i,j-1] Ω[i,j+1] Ω[i-1,j] Ω[i+1,j]]
    
    if     a == [0 1 2 2]
        n = [-1 0]
    elseif a == [1 0 2 2]
        n = [1 0]
    elseif a == [2 2 1 0]
        n = [0 -1]
    elseif a == [2 2 0 1]
        n = [0 1]
    elseif a == [0 2 0 2]
        n = [-1 1]
    elseif a == [2 0 0 2]
        n = [1 1]
    elseif a == [0 2 2 0]
        n = [-1 -1]
    elseif a == [2 0 2 0]
        n = [1 -1]
    else n = [0 0]
    end
end

#nodos_neumann() modifica los puntos imaginarios de la frontera neumann.
#de esta forma, podemos usar la ec de discretización para u(i,j,n+1)
#de acuerdo a la normal n
#el nodo Ω_ij

function nodos_neumann(U,i,j,n,R)   ####Coef de Reflexion R
    if n[1] == -1
        U[i,j-1] = R*U[i,j+1]
    elseif n[1] == 1
        U[i,j+1] = R*U[i,j-1]
    end
    
    if n[2] == -1
        U[i+1,j] = R*U[i-1,j]
    elseif n[2] == 1
        U[i-1,j] = R*U[i+1,j]
    end
    return U
end


#Función tipo heaviside dentro de una elipse
function q(x,y,k,l)
    if (x/k)^2 + (y/l)^2 <= 1/4
        q = 1
    else
        q = 0
    end
end

#el nodo de la malla de Q más cercano a la posición  de un agente 
#es donde se localizara la fuente "discreta" de contaminantes.
function Q_gen(X,Y,x)
    N,M = length(X),length(Y)
    k,l = X[2]-X[1], Y[2]-Y[1]
    Q = zeros(M,N)
    for i in 1:M
        for j in 1:N
            Q[i,j] = q(X[j]-x[1],Y[i]-x[2],k,l)
        end
    end
    return Q
end

# Qα = [Q_gen(xt1),... , Q_gen(xtn) ] calcula la matriz Q de un agente
# en tiempos t_1,...,t_n

function Qα(X, Y, Xα)   #Xα =[x(t1),x(t2),...,x(tn)]
    Qα = []
    for i in 1: length(Xα)
        push!( Qα, Q_gen(X,Y,Xα[i]) )
    end
    return Qα
end

# QT contiene la distribución de fuentes de contaminantes
# en tiempos t_1,...,t_n

function QT(X, Y, XT)   #XT = [Xα1,... XαN]
    QT = Qα(X,Y, XT[1])
    if length(XT) > 1
        for i in 2:length(XT)
            QT = QT + Qα(X,Y,XT[i])
        end
    end
    return QT
end

function dif_ij(U,Q,i,j,D,h,k,l)
    u_ij = D*h*((U[i,j+1] - 2*U[i,j] + U[i,j-1])/(k^2) + (U[i+1,j] - 2*U[i,j] + U[i-1,j])/(l^2) ) + h*Q[i,j] + U[i,j]
end

function paso_dif(u,Ω,X,Y,D,h,Q,R)
    u = rodear(u,0)
    Ω = rodear(Ω,0)
    Q = rodear(Q,0)
    M, N = size(u)
    k,l = Float32.(X[2]-X[1]), Float32.(Y[2]-Y[1])
    u2 = Float32.(zeros(M,N))
    
    if h <= k^2*l^2/(2*D*(k^2 +l^2))
    
    for i in 1:M
        for j in 1:N
            if Ω[i,j] == 1
                u2[i,j] = dif_ij(u,Q,i,j,D,h,k,l)
            elseif Ω[i,j] == 2
                n = id_normal(Ω,i,j)
                u = nodos_neumann(u,i,j,n,R)
                u2[i,j] = dif_ij(u,Q,i,j,D,h,k,l)
                
            elseif Ω[i,j] == 0
                u2[i,j] = 0
            end
        end
    end
        
    else println("no se satisfacen condiciones de estabilidad")
        @show h,k^2*l^2/(2*D*(k^2 +l^2))
    end
    
    return u2[2:M-1, 2:N-1]
end

#Calcula la difusión para n pasos de tiempo

function difusion_mean(u,Ω,X,Y,D,h,XT,R,n)
    
    #QTT = QT(X, Y, XT)
    U = Float32.(u)
    UU = Float32.([])
    M = []
    
    for i in 1:n
        
        XTi = [XT[j][i] for j in 1:length(XT)]
        QTT = Float32.(sum(Qα(X, Y, XTi)))
        U = paso_dif(U,Ω,X,Y,D,h,QTT,R)
        push!(M,mean(U))
    end
    return U,M
end

#Calcula la difusión para n pasos de tiempo

function difusion(u,Ω,X,Y,D,h,XT,R,n)
    
    #QTT = QT(X, Y, XT)
    U = u
    UU = []
    M = []
    
    for i in 1:n
        
        XTi = [XT[j][i] for j in 1:length(XT)]
        QTT = Float32.(sum(Qα(X, Y, XTi)))
        U = paso_dif(U,Ω,X,Y,D,h,QTT,R)
        #if mod(i,200) == 1
         #   push!(UU,U)
        #end
        push!(UU,U)
    end
    return UU
end

#Calcula el promendio de una lista de valores en el arreglo A
function multiply(a)
    m = 1
    for i in 1:length(a)
        m = m*a[i]
    end
    return m
end

function mean(A)
    mean = sum(A) /multiply(size(A))
end

#contiene el valor de u en el nodo i,j de las n iteraciones
function u_punto_xy(UU,j,i)
    n = length(UU)
    return [(UU[k])[j,i] for k in 1:n]
end

#Distribución del valor U al tiempo t
function distribucion_t(UU,t)
    M, N = size(UU[t])
    D_x = sum([UU[t][j,:] for j in 1:M])
    #D_x = D_x/(sum(D_x))
    D_y = sum([UU[t][:,i] for i in 1:N])/N
    #D_y = D_y/(sum(D_y))
    
    return D_x,D_y
end

function norm(v)
    l = length(v)
    S = 0
    for i in 1:l
        S += v[i]^2
    end
    return sqrt(S)
end

function histoU(U,n,min,max)
    M, N = size(U)
    h = zeros(n)
    l = (max-min)/n
    for i in 1:M
        for j in 1:N
            w = Int(floor(U[i,j]/l))+1
            h[w] = h[w] + 1
        end
    end
    
    return h/norm(h)
end
