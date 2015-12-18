for i=1:Zsize

% Full Cell Integrated Intensity
    FCII{i}=ydatacrpd{i}(lob:upb,1:XSize(i));
    x(i,7)=sum(sum(FCII{i}));
    
    % SPOT around centroid
    if SX(i)==1 || SX(i)==2
        
        % shift x when it is too close to the left border
        if PX(i)==0 || PX(i)==1  
        x(i,2)=2;
        end
        
        % shift x when it is too close to the right border
        if PX(i)==XSize(i) || PX(i)==XSize(i)-1 
        x(i,2)=XSize(i)-2;
        end
        
        
    for j=[-1 0 1] 
        if round(x(i,5))==2 || round(x(i,5))==1
        for k=[-1 0 1]
    II{j+2,k+2,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+j);
        end
        end
        if round(x(i,5))==3 || round(x(i,5))==4
        for k=[-2 -1 0 1 2]
    II{j+2,k+3,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+j);
        end
        end
    end
    end
    
    if round(x(i,3))==3 || round(x(i,3))==4
        
        % shift x when it is too close to the left border
        if round(x(i,2))==0 || round(x(i,2))==1 || round(x(i,2))==2 
        x(i,2)=3;
        end
        
        % shift x when it is too close to the right border
        if round(x(i,2))==XSize(i) || round(x(i,2))==XSize(i)-1 || round(x(i,2))==XSize(i)-2 
        x(i,2)=XSize(i)-3;
        end
        
        
    for j=[-2 -1 0 1 2] 
        if round(x(i,5))==2 || round(x(i,5))==1
        for k=[-1 0 1]
    II{j+3,k+2,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+j);
        end
        end
        if round(x(i,5))==3 || round(x(i,5))==4
        for k=[-2 -1 0 1 2]
    II{j+3,k+3,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+j);
        end
        end
    end
    end
    
        SII=[II{:,:,i}];
        x(i,6)=sum(SII);

%---------- Second Spot ---------------------------------------------------
    if round(xR1(i,3))==2 || round(xR1(i,3))==1
        
        % shift x when it is too close to the left border
        if round(xR1(i,2))==0 || round(xR1(i,2))==1  
        xR1(i,2)=2;
        end
        
        % shift x when it is too close to the right border
        if round(xR1(i,2))==XSize(i) || round(xR1(i,2))==XSize(i)-1 
        xR1(i,2)=XSize(i)-2;
        end
        
        
    for j=[-1 0 1] 
        if round(xR1(i,5))==2 || round(xR1(i,5))==1
        for k=[-1 0 1]
    II1{j+2,k+2,i}=ydatacrpdR1{i}(round(xR1(i,4))+k,round(xR1(i,2))+j);
        end
        end
        if round(xR1(i,5))==3 || round(xR1(i,5))==4
        for k=[-2 -1 0 1 2]
    II1{j+2,k+3,i}=ydatacrpdR1{i}(round(xR1(i,4))+k,round(xR1(i,2))+j);
        end
        end
    end
    end
    
    if round(xR1(i,3))==3 || round(xR1(i,3))==4
        
        % shift x when it is too close to the left border
        if round(xR1(i,2))==0 || round(xR1(i,2))==1 || round(xR1(i,2))==2 
        xR1(i,2)=3;
        end
        
        % shift x when it is too close to the right border
        if round(xR1(i,2))==XSize(i) || round(xR1(i,2))==XSize(i)-1 || round(xR1(i,2))==XSize(i)-2 
        xR1(i,2)=XSize(i)-3;
        end
        
        
    for j=[-2 -1 0 1 2] 
        if round(xR1(i,5))==2 || round(xR1(i,5))==1
        for k=[-1 0 1]
    II1{j+3,k+2,i}=ydatacrpdR1{i}(round(xR1(i,4))+k,round(xR1(i,2))+j);
        end
        end
        if round(xR1(i,5))==3 || round(xR1(i,5))==4
        for k=[-2 -1 0 1 2]
    II1{j+3,k+3,i}=ydatacrpdR1{i}(round(xR1(i,4))+k,round(xR1(i,2))+j);
        end
        end
    end
    end
    
        SII1=[II1{:,:,i}];
        xR1(i,6)=sum(SII1);
        
        
%---------- Third Spot ---------------------------------------------------
    
    if round(xR2(i,3))==2 || round(xR2(i,3))==1
        
        % shift x when it is too close to the left border
        if round(xR2(i,2))==0 || round(xR2(i,2))==1  
        xR2(i,2)=2;
        end
        
        % shift x when it is too close to the right border
        if round(xR2(i,2))==XSize(i) || round(xR2(i,2))==XSize(i)-1 
        xR2(i,2)=XSize(i)-2;
        end
        
        
    for j=[-1 0 1] 
        if round(xR2(i,5))==2 || round(xR2(i,5))==1
        for k=[-1 0 1]
    II2{j+2,k+2,i}=ydatacrpdR2{i}(round(xR2(i,4))+k,round(xR2(i,2))+j);
        end
        end
        if round(xR2(i,5))==3 || round(xR2(i,5))==4
        for k=[-2 -1 0 1 2]
    II2{j+2,k+3,i}=ydatacrpdR2{i}(round(xR2(i,4))+k,round(xR2(i,2))+j);
        end
        end
    end
    end
    
    
    if round(xR2(i,3))==3 || round(xR2(i,3))==4
        
        % shift x when it is too close to the left border
        if round(xR2(i,2))==0 || round(xR2(i,2))==1 || round(xR2(i,2))==2 
        xR2(i,2)=3;
        end
        
        % shift x when it is too close to the right border
        if round(xR2(i,2))==XSize(i) || round(xR2(i,2))==XSize(i)-1 || round(xR2(i,2))==XSize(i)-2 
        xR2(i,2)=XSize(i)-3;
        end
        
        
    for j=[-2 -1 0 1 2] 
        if round(xR2(i,5))==2 || round(xR2(i,5))==1
        for k=[-1 0 1]
    II2{j+3,k+2,i}=ydatacrpdR2{i}(round(xR2(i,4))+k,round(xR2(i,2))+j);
        end
        end
        if round(xR2(i,5))==3 || round(xR2(i,5))==4
        for k=[-2 -1 0 1 2]
    II2{j+3,k+3,i}=ydatacrpdR2{i}(round(xR2(i,4))+k,round(xR2(i,2))+j);
        end
        end
    end
    end
    
        SII2=[II2{:,:,i}];
        xR2(i,6)=sum(SII2);
        
%---------- Fourth Spot ---------------------------------------------------

    if round(xR3(i,3))==2 || round(xR3(i,3))==1

        % shift x when it is too close to the left border
        if round(xR3(i,2))==0 || round(xR3(i,2))==1  
        xR3(i,2)=2;
        end
        
        % shift x when it is too close to the right border
        if round(xR3(i,2))==XSize(i) || round(xR3(i,2))==XSize(i)-1 
        xR3(i,2)=XSize(i)-2;
        end
        
        
    for j=[-1 0 1] 
        if round(xR3(i,5))==2 || round(xR3(i,5))==1
        for k=[-1 0 1]
    II3{j+2,k+2,i}=ydatacrpdR3{i}(round(xR3(i,4))+k,round(xR3(i,2))+j);
        end
        end
        if round(xR3(i,5))==3 || round(xR3(i,5))==4
        for k=[-2 -1 0 1 2]
    II3{j+2,k+3,i}=ydatacrpdR3{i}(round(xR3(i,4))+k,round(xR3(i,2))+j);
        end
        end
    end
    end
    
    if round(xR3(i,3))==3 || round(xR3(i,3))==4
    
        % shift x when it is too close to the left border
        if round(xR3(i,2))==0 || round(xR3(i,2))==1 || round(xR3(i,2))==2 
        xR3(i,2)=3;
        end
        
        % shift x when it is too close to the right border
        if round(xR3(i,2))==XSize(i) || round(xR3(i,2))==XSize(i)-1 || round(xR3(i,2))==XSize(i)-2 
        xR3(i,2)=XSize(i)-3;
        end
    
    for j=[-2 -1 0 1 2] 
        if round(xR3(i,5))==2 || round(xR3(i,5))==1
        for k=[-1 0 1]
    II3{j+3,k+2,i}=ydatacrpdR3{i}(round(xR3(i,4))+k,round(xR3(i,2))+j);
        end
        end
        if round(xR3(i,5))==3 || round(xR3(i,5))==4
        for k=[-2 -1 0 1 2]
    II3{j+3,k+3,i}=ydatacrpdR3{i}(round(xR3(i,4))+k,round(xR3(i,2))+j);
        end
        end
    end
    end
    
        SII3=[II3{:,:,i}];
        xR3(i,6)=sum(SII3(:));
end