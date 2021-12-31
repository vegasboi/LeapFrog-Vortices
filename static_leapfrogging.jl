# Calculates hypotanuse with 2 float arguments
function hyp(a,b)
        return(a*a+b*b)^0.5
       end
# Calculates the distance between 2 points (x1,y1) and (x2,y2)
function dist(x1,y1,x2,y2)
        return hyp(x1-x2,y1-y2)
       end
# Returns induced velocity vector of point (x1,y1) from positioning of point (x2,y2)
function vector(x1,y1,x2,y2)
        return (x1-x2)/(2*pi*dist(x1,y1,x2,y2)^2),(y1-y2)/(2*pi*dist(x1,y1,x2,y2)^2)
       end
using LinearAlgebra
#= Returns induced velocity vector of points (x1,y1) (x2,y2) (x3,y3) (x4,y4) given 
   their respective positionings =#
 function velocity(x1,y1,x2,y2,x3,y3,x4,y4)
        gamma_top = [0,0,-1]     #CCW rotation of the top vortices
        gamma_bottom = [0,0,1]   #CW rotation of the bottom vortices

        vor1_vec1_x,vor1_vec1_y = vector(x1,y1,x2,y2) #Normalized 2-d vector between (x1,y1) and (x2,y2)
        vor1_vec1 = [vor1_vec1_x,vor1_vec1_y,0] #Normalized 3-d vector between (x1,y1) and (x2,y2)
        vor1_vec1 = cross(vor1_vec1,gamma_top) #Induced velocity vectory on (x1,y1) from (x2,y2)

        vor1_vec2_x,vor1_vec2_y = vector(x1,y1,x3,y3) #Normalized 2-d vector between (x1,y1) and (x3,y3)
        vor1_vec2 = [vor1_vec2_x,vor1_vec2_y,0]
        vor1_vec2 = cross(vor1_vec2,gamma_top)

        vor1_vec3_x,vor1_vec3_y = vector(x1,y1,x4,y4) #Normalized 2-d vector between (x1,y1) and (x4,y4)
        vor1_vec3 = [vor1_vec3_x,vor1_vec3_y,0]
        vor1_vec3 = cross(vor1_vec3,gamma_bottom)

        vor1_vec = vor1_vec1 + vor1_vec2 + vor1_vec3 #Resultant of all the induced velocities on point (x1,y1)
        vor2_vec = [vor1_vec[1],-1*vor1_vec[2],0] #Resultant of all the induced velocities on point (x2,y2)
        
        #Calculating the induced velocity on point (x3,y3) and (x4,y4)
        vor3_vec1_x,vor3_vec1_y = vector(x3,y3,x1,y1) 
        vor3_vec1 = [vor3_vec1_x,vor3_vec1_y,0]
        vor3_vec1 = cross(vor3_vec1,gamma_bottom)

        vor3_vec2_x,vor3_vec2_y = vector(x3,y3,x2,y2)
        vor3_vec2 = [vor3_vec2_x,vor3_vec2_y,0]
        vor3_vec2 = cross(vor3_vec2,gamma_top)

        vor3_vec3_x,vor3_vec3_y = vector(x3,y3,x4,y4)
        vor3_vec3 = [vor3_vec3_x,vor3_vec3_y,0]
        vor3_vec3 = cross(vor3_vec3,gamma_bottom)

        vor3_vec = vor3_vec1 + vor3_vec2 + vor3_vec3
        vor4_vec = [vor3_vec[1],-1*vor3_vec[2],0]

        return vor1_vec[1],vor1_vec[2],vor2_vec[1],vor2_vec[2],vor3_vec[1],vor3_vec[2],vor4_vec[1],vor4_vec[2]
       end
 x1_pos = 0.00  #Bottom left vortice initial x position
 y1_pos = -0.500 #Bottom left vortice initial y position
 x2_pos = 0.00 #Top left vortice initial x position
 y2_pos = 0.500 #Top left vortice initial y position
 x3_pos = 1.00 #Top right vortice initial x position
 y3_pos = 0.500 #Top right vortice initial y position
 x4_pos = 1.00 #Bottom right vortice initial x position
 y4_pos = -0.500 #Bottom right vortice initial y position
 #Creating vectors to store x and y positions of each vortice over time
 x1 = [x1_pos]
 y1 = [y1_pos]
 x2 = [x2_pos]
 y2 = [y2_pos]
 x3 = [x3_pos]
 y3 = [y3_pos]
 x4 = [x4_pos]
 y4 = [y4_pos] 
#Tracking vortice trajectories over 40 seconds with a time step of 0.01 seconds
for i in 0:0.01:40
        #Using discrete approximations to track trajectories with updated velocity vectors and positions
        x1_vel,y1_vel,x2_vel,y2_vel,x3_vel,y3_vel,x4_vel,y4_vel = velocity(x1_pos,y1_pos,x2_pos,y2_pos,x3_pos,y3_pos,x4_pos,y4_pos)
        global x1_pos = x1_pos + x1_vel*0.01
        global y1_pos = y1_pos + y1_vel*0.01
        global x2_pos = x2_pos + x2_vel*0.01
        global y2_pos = y2_pos + y2_vel*0.01
        global x3_pos = x3_pos + x3_vel*0.01
        global y3_pos = y3_pos + y3_vel*0.01
        global x4_pos = x4_pos + x4_vel*0.01
        global y4_pos = y4_pos + y4_vel*0.01
        
        push!(x1,x1_pos)
        push!(y1,y1_pos)
        push!(x2,x2_pos)
        push!(y2,y2_pos)
        push!(x3,x3_pos)
        push!(y3,y3_pos)
        push!(x4,x4_pos)
        push!(y4,y4_pos)    
        println(x1_pos,y1_pos,x2_pos,y2_pos)
       end
using Plots
p = plot(x1,y1)
plot!(p,x2,y2)
plot!(p,x3,y3)
plot!(p,x4,y4)
plot!(xlims = (0,11), ylims = (-2,2),xlabel = "Horizontal Displacement (cm)", ylabel = "Vertical Displacement (cm)", title = "Vortex Trajectory")
