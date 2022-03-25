#
# Array iteration
#

# Matriz de 2 filas y 3 columnas
#
m = rand(2,3)

#=
 Recorrer la matriz por filas y columnas
 al poner el indice j en el loop interno primero recorremos las filas
=#
for i in 1:size(m,1)
    for j in 1:size(m,2)
        @info "$(i) - $(j) = $(m[i,j])"
    end
end
m

#=
 Recorrer la matriz por posición en la memoria
 las matrices están guardadas por columnas
=#
for i in eachindex(m)
    @info "$(i) = $(m[i])"
end

#=
 Lo mismo que antes pero usando enumerate 
 esta sería la forma más optimizada de iterar por la matriz
=#
for (i, value ) in enumerate(m)
    @info "$(i) = $(value)"
end

