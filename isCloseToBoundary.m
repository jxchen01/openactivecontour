function flag=isCloseToBoundary(P,dimx,dimy)

bufferSize = 4;
flag = false;

if( any(P([1,end],1) > dimx - bufferSize) || any(P([1,end],1)<=bufferSize) )
    flag = true;
    return
end

if( any(P([1,end],2) > dimy - bufferSize) || any(P([1,end],2)<=bufferSize) )
    flag = true;
    return
end
