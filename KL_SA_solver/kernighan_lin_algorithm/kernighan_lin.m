function [P1, P2, cutValue] = kernighan_lin(W, P1, P2)
%%   Kernighan-Lin algoritam za pronalazenje minimalnog cutset-a

%   Problem: particionisanje el.kola na dva dela tako da se preseca
%   najmanji broj veza izmedju logickih kola (problem rutiranja)
%   W - matrica susednosti, sa tezinama svake grane
%   P1, P2 - inicijalne particije
%   Povratna vrednost funkcije su particije cvorova P1, P2 i velicina
%   cut size

    n = size(W,1); 
    improved = true;
    it = 0;it1 = 0;
    while improved %kada se nadje prvi slucaj gde se rezultat pogorsao, prekida se alg.
    it = it+1;    
    
        improved = false;
        locked = false(1,n);%na pocetku su svi otkljucani 
        gains = zeros(1, floor(n/2));
        pairList = zeros(floor(n/2), 2);
        
        
        D = computeD(W, P1, P2);
        
        % 
         for k = 1:floor(n/2)
            
            maxGain = -inf;
            bestPair = [-1 -1];

           %postavke a i b su inicijalni - fiksirani u prvom prolazu
           if k==1 
            maxGain = D(P1(1)) + D(P2(1)) - 2*W(1,2);
            bestPair = [1 2];
           
           else

            % Pronadji otključan par čvorova sa najvecim gain-om 
            for ai = 1:length(P1)
                a = P1(ai);
                if locked(a), continue; end
                for bi = 1:length(P2)
                    b = P2(bi);
                    if locked(b), continue; end
                    g = D(a) + D(b) - 2*W(a,b); %računanje gain-a
                    if g > maxGain
                        maxGain = g;
                        bestPair = [a b];
                    end
                end
            end
            
           end
            % Zaključaj pronadjen par
            a = bestPair(1); b = bestPair(2);
            locked([a b]) = true;
            gains(k) = maxGain;
            pairList(k,:) = [a b];
            
            % D vrednosti otkljucanih parova
            for v = 1:n
                if ~locked(v)
                    if ismember(v,P1)
                        D(v) = D(v) + 2*(W(v,b) - W(v,a));
                    else
                        D(v) = D(v) + 2*(W(v,a) - W(v,b));
                    end
                end
            end
        end
        
        % Izracunaj trenutni cut size i cut size nakon zamene
        currentCut = computeCut(W, P1, P2);
        newP1 = P1; newP2 = P2;
        %provera drugih mogucih kombinacija
        for k = 1:floor(n/2)
            a = pairList(k,1);
            b = pairList(k,2);
            %swap
            newP1(newP1==a) = [];
            newP2(newP2==b) = [];
            newP1 = [newP1 b];
            newP2 = [newP2 a];
            
            newCut = computeCut(W, newP1, newP2);
       %Da li se rezultat poboljsao?
            if newCut < currentCut
                % Jeste
            %disp('_____________________________________')    
            disp('nova particija P1');disp(newP1)
            disp('nova particija P2');disp(newP2)
            fprintf('Novi cut size %d \n', newCut);
            disp('_____________________________________')    

                currentCut = newCut;
                improved = true;
                P1 = newP1;
                P2 = newP2;
            else
                % Nije, vrati prethodno
                newP1(newP1==b) = [];
                newP2(newP2==a) = [];
                newP1 = [newP1 a];
                newP2 = [newP2 b];
            end
        end
    
    
    end
    
    % Compute final cut size
    cutValue = computeCut(W, P1, P2);
    
end
%--------------------------------------------------------------------------
%   Funkcija koja racuna razliku D : external_cost - internal_cost za dati
%   cvor, 
%   ext -> suma tezina grana koje povezuju cvor x i njemu susedne u drugoj
%   particiji; 
%   int -> suma tezina grana koje povezuju cvor x i njemu susedne u
%   particiji kojoj pripada x


function D = computeD(W, A, B)
    n = size(W,1);
    D = zeros(1,n);
    for v = 1:n
        ext = 0; int = 0;

        if ismember(v,A) % ako cvor pripada A
            ext = sum(W(v,B));
            int = sum(W(v,A));
        else                % ako cvor pripada B
            ext = sum(W(v,A));
            int = sum(W(v,B));
        end

        D(v) = ext - int;
    end
end
%--------------------------------------------------------------------------
%   Funkcija koja racuna cutsize: broj grana koje preseca linija koja deli
%   graf na dve particije

function cutVal = computeCut(W, A, B)
   
    cutVal = nnz(W(A,B));  % nnz = number of nonzero elements
end
