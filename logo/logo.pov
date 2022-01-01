global_settings{assumed_gamma 1.0}

#declare Lng=20;
#declare Lat=20;
#declare Tilt=0;
#declare Pers=0.1;
#declare Zoom=0.65;
#declare LookAt=<0,0,0>;

#macro SCS(lng,lat) <cos(radians(lat))*cos(radians(lng)),cos(radians(lat))*sin(radians(lng)),sin(radians(lat))> #end
#declare AspectRatio=image_width/image_height;
#declare Z=SCS(Lng,Lat);
#declare X=vaxis_rotate(<-sin(radians(Lng)),cos(radians(Lng)),0>,Z,Tilt);
#declare Y=vcross(Z,X);
#if(Pers)
    #declare Loc=LookAt+SCS(Lng,Lat)/(Zoom*Pers);
    camera{
        perspective
        location Loc
        right -2*X*sqrt(AspectRatio)/Zoom
        up 2*Y/(sqrt(AspectRatio)*Zoom)
        direction Z/(Zoom*Pers)
        sky Y
        look_at LookAt
    }
    light_source{
        Loc
        color rgb<1,1,1>*2
    }
#else
    #declare Loc=SCS(Lng,Lat);
    camera{
        orthographic
        location Loc*100
        right -2*X*sqrt(AspectRatio)/Zoom
        up 2*Y/(sqrt(AspectRatio)*Zoom)
        sky Y
        look_at LookAt
    }
    light_source{
        SCS(Lng,Lat)
        color rgb<1,1,1>
        parallel
        point_at 0
    }
#end
background{rgb<1,1,1>}

// cylinder{<0,0,0>,<1,0,0>,0.1 pigment{rgb<1,0,0>}}
// cylinder{<0,0,0>,<0,1,0>,0.1 pigment{rgb<0,1,0>}}
// cylinder{<0,0,0>,<0,0,1>,0.1 pigment{rgb<0,0,1>}}

#macro Rotate(w,i)
    #local j=i;
    #local wj=w;
    #while(j)
        #local wj=<wj.z,wj.x,wj.y>;
        #local j=j-1;
    #end
    wj
#end

#declare i=-1;
#while(i<2)
    #declare j=-1;
    #while(j<2)
        #declare k=-1;
        #while(k<2)
            sphere{<i,j,k>,0.01 pigment{rgb 0.05}}
            #declare k=k+2;
        #end
        #declare l=0;
        #while (l<3)
            cylinder{Rotate(<i,j,1>,l),Rotate(<i,j,-1>,l),0.01 pigment{rgb 0.05}}
            #declare l=l+1;
        #end
        #declare j=j+2;
    #end
    #declare i=i+2;
#end

#macro F(X,Y)
    <X, Y, X*X+Y*Y-1>
#end

#declare paraboloid = union{
#declare N = 50;
#declare i=0;
#while(i<N)
    #declare j=0;
    #while(j<N)
        #declare d = 2/N;
        #declare x_ = -1+i*d;
        #declare y_ = -1+j*d;
        #declare p0 = F(x_+d/2,y_+d/2);
        #declare p1 = F(x_,y_);
        #declare p2 = F(x_+d,y_);
        #declare p3 = F(x_+d,y_+d);
        #declare p4 = F(x_,y_+d);
        triangle{p1,p2,p0}
        triangle{p2,p3,p0}
        triangle{p3,p4,p0}
        triangle{p4,p1,p0}
        #declare j=j+1;
    #end
    #declare i=i+1;
#end
}

// object{paraboloid pigment{image_map{jpeg "/home/hyrodium/Git/ElasticSurfaceEmbedding/img/me.jpg" map_type 0}}}

object{paraboloid pigment{image_map{png "ura.png" map_type 0} scale 2}}
object{paraboloid pigment{image_map{png "omote.png" map_type 0} scale 2} translate <0, 0, -0.001>}
