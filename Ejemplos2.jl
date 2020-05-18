include("Difusion_funciones.jl")

using Plots

N,M = 101,101
X,Y = linspace(-10,10,N), linspace(-10,10,M)
h,D = 0.01,0.1

u = rand(M,N)
U = u
Ω = rodear(ones(M-2,N-2),2)
Q = zeros(M,N)
UU = []
n=500

k,l = X[2]-X[1], Y[2]-Y[1]
h <= k^2*l^2/(2*D*(k^2 +l^2)), h, k^2*l^2/(2*D*(k^2 +l^2)) #condicion de estabilidad

for i in 1:n
    U = paso_dif(U,Ω,X,Y,D,h,Q,1.0)
    push!(UU,U)
end

t = 500
surface(X,Y,UU[t],zlim=[0.2,0.8],leg=false, title= "t=$t", aspectratio=:equal,color=:magma,levels=25)

plot(linspace(0,n*h,n),[mean(UU[i]) for i in 1:n], xlabel= "t",ylabel="Valor esperado de u",grid=false)

#Discretización

N,M =201,101
X = linspace(0,200,N)
Y = linspace(-40,40,M)
Q = zeros(M,N)
h = 1.0
D = 0.15


#condiciones iniciales
u = zeros(M,N)
U = u
UU = []

#iteraciones
n = 6000


#estabilidad
k,l = X[2]-X[1], Y[2]-Y[1]
h <= k^2*l^2/(2*D*(k^2 +l^2)), h, k^2*l^2/(2*D*(k^2 +l^2))

#Definición del espacio Ω

Ω = ones(M-2,N-2)
Ω = rodear(Ω,0)
Block2 = rodear(zeros(4,10),2)
Block0 = zeros(6,14)

for i in 1:30
    a,b = rand(1:M-5),rand(1:N-12)
    Ω[a:a+5,b:b+11] = Block2
end

for i in 1:20
    a,b = rand(1:M-5),rand(1:N-13)
    Ω[a:a+5,b:b+13] = Block0
end

contour(X,Y,Ω)

#generar trayectoria XT

x(t) = [mod(0.8*t,200),0]
X1 = []
for i in 1:n+1
    push!(X1, x(i-1))
end
XT = [X1];

UU = difusion(u,Ω,X,Y,D,h,XT,0.8,n);

contour(X,Y,Ω, color=:thermal,xlabel="x", ylabel="y")
contour!(X,Y,UU[400],leg=false,levels=150, color=:magma,zlim=[0,1],grid=false)

plot([mean(UU[i]) for i in 1:n],xlabel="tiempo", ylabel="<u>", w=2,grid=false, color=:grey)

plot(Y,distribucion_t(UU,200)[2], w=2,xlabel ="y",label="t = 200",ylabel="Distribución de contaminantes")
plot!(Y,distribucion_t(UU,1000)[2],w=2, xlabel ="y",label="t = 1000",ylabel="Distribución de contaminantes")
plot!(Y,distribucion_t(UU,4000)[2], w=2,xlabel ="y",label="t = 4000",ylabel="Distribución de contaminantes",grid=false)


