/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof.: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 ******************************************************************************/

/* IMPLEMENTE AS FUNCOES INDICADAS
 * DIGITE SEU RA: 145980
 * SUBMETA SOMENTE ESTE ARQUIVO */

#include <iostream>
#include <float.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "lpdtspalgs.h"

bool naive(const LpdTspInstance &l, LpdTspSolution  &s, int tl);

bool constrHeur(const LpdTspInstance &l, LpdTspSolution  &s, int tl)
/* Implemente esta função, entretanto, não altere sua assinatura */
{
	// solucao inspirada em algoritmo guloso
	
	double capacity = l.capacity;
	double available = capacity;
	DNode deposito = l.depot;
	
	vector<int> transportado (l.k);
	
	int i;
	for(i = 0; i < l.k; i++){
		// 0 - nao transportado
		// 1 - em posse do caixeiro
		// 2 - entregue
		transportado[i] = 0;
	}
	
	s.tour.clear();
	s.cost = 0.0;
	
	// insere o deposito no tour
	s.tour.push_back(deposito);
	
	// procura o vizinho mais proximo do deposito que
	// contenha um item a ser recolhido
	DNode vizinhoMaisProx = deposito;
	double minWeight = DBL_MAX;
	for(OutArcIt e(l.g, deposito); e != INVALID; ++e){
		if(l.weight[e] < minWeight && l.s[l.g.target(e)] > 0){
			minWeight = l.weight[e];
			vizinhoMaisProx = l.g.target(e);
		}
	}
	
	// insere o vizinho mais proximo no tour
	DNode prox = vizinhoMaisProx;
	s.tour.push_back(prox);
	s.cost += minWeight;
	
	for(;;){ // a interrupcao deste for eh por break durante a execucao
		
		// se ha itens a serem entregues e que estejam com o caixeiro
		// entao serao entregues
		if( l.t[prox] > 0 ){
			int idItem = l.t[prox];
			if(transportado[idItem - 1] == 1 ){
				transportado[idItem - 1] = 2;
				available = available + l.items[idItem-1].w;
			}
		}
		
		// se ha itens a serem recolhidos e que caibam na capacidade
		// do caixeiro, entao serao recolhidos
		if( l.s[prox] > 0 ){
			int idItem = l.s[prox];
			if(transportado[idItem - 1] == 0 && ( (available - l.items[idItem-1].w) >= 0 ) ){
				transportado[idItem - 1] = 1;
				available = available - l.items[idItem-1].w;
			}
		}
		
		/*
		for(i = 0; i < l.k; i++){
			printf("transportado[%d] = %d\n", i, transportado[i]);
		}
		printf("-----------------\n");
		*/
		
		// verifica se ja entregou todos os itens
		// caso positivo, eh necessario interromper o for
		bool continua = false;
		for(i = 0; i < l.k; i++){
			if(transportado[i] != 2){
				continua = true;
			}
		}
		
		if(!continua){
			break;
		}
		
		// busca um proximo vertice a ser visitado
		// que seja promissor no sentido de economia de custos
		double minWeight1 = DBL_MAX; // menor distancia para entrega
		double minWeight2 = DBL_MAX; // menor distancia para recolhimento viavel
		DNode aux1 = deposito;
		DNode aux2 = deposito;
		Arc e;
		for(i = 0; i < l.k; i++){
			if(transportado[i] == 1){
				e = findArc(l.g, prox, l.items[i].t);
				if(l.weight[e] < minWeight1){
					aux1 = l.items[i].t;
					minWeight1 = l.weight[e];
				}
			}

			if(transportado[i] == 0){
				e = findArc(l.g, prox, l.items[i].s);
				if(l.weight[e] < minWeight2 && ( (available - l.items[i].w) >= 0 )){
					aux2 = l.items[i].s;
					minWeight2 = l.weight[e];
				}
			}
		}

		// faz o movimento que parece ser mais economico
		if(minWeight1 != DBL_MAX || minWeight2 != DBL_MAX){ // verifica se achou arco
			if(minWeight1 < minWeight2){
				prox = aux1;
				s.cost += minWeight1;
				s.tour.push_back(prox);
			}else{
				prox = aux2;
				s.cost += minWeight2;
				s.tour.push_back(prox);
			}
		}else{ // se nao achou vai simplesmente para o vertice mais proximo
			minWeight = DBL_MAX;
			vizinhoMaisProx = deposito;
			for(OutArcIt e(l.g, prox); e != INVALID; ++e){
				if(l.weight[e] < minWeight && l.s[l.g.target(e)] > 0){
					minWeight = l.weight[e];
					vizinhoMaisProx = l.g.target(e);
				}
			}
			prox = vizinhoMaisProx;
			s.cost += minWeight;
			s.tour.push_back(prox);
		}
	}
	
	// por fim, insere o arco que liga o ultimo vertice do tour ao deposito
	Arc e = findArc(l.g, prox, deposito);
	s.cost += l.weight[e];
	s.lowerBound = s.cost;
	
	return false;
}
//------------------------------------------------------------------------------
bool metaHeur(const LpdTspInstance &l, LpdTspSolution  &s, int tl)
/* Implemente esta função, entretanto, não altere sua assinatura */
{
	// BUSCA TABU
	
	clock_t begin = clock();
	
	// utiliza a heuristica construtiva como solucao inicial
	constrHeur(l, s, tl);
	
	clock_t posCons = clock();
	if( (double(posCons - begin) / CLOCKS_PER_SEC) > tl ){
		return false;
	}
	
	int x = 1;
	
	clock_t loopClock = clock();
	do{		
		int i, j;
		int max = s.tour.size();
		for(i = 0; i < max-x; i++){
			vector<DNode> tourEval(max);
			for(j = 0; j < max; j++){
				tourEval[j] = s.tour[j];
			}
			
			if( i+x < (int)tourEval.size() ){
				DNode temp = tourEval[i];
				tourEval[i] = tourEval[i+x];
				tourEval[i+x] = temp;
				
				bool valido = true;			
				
				// aqui foram utilizadas linhas de codigo da validacao
				// de solucao por Mauro H. Mulati
				vector<bool> picked(l.k, false);         // 0..k-1 to represent items 1..k. Initialized as none of the items is picked
				vector<bool> delive(l.k, false);         // 0..k-1 to represent items 1..k. Initialized as none of the items is delivered
				for(int i = 1; i < (int)tourEval.size(); i++){  // Do not look at the depot (first node)
				  if(l.s[ tourEval[i] ] > 0){           // If in a DNode a item starts
					 if(!picked[ l.s[ tourEval[i] ] - 1 ] && !delive[ l.s[ tourEval[i] ] - 1 ]){
						picked[ l.s[ tourEval[i] ] - 1 ] = true;
					 }
					 else{
						valido = false;
					 }
				  }
				  else if(l.t[ tourEval[i] ] > 0){      // If in a DNode a item terminates
					 if(picked[ l.t[ tourEval[i] ] - 1 ] && !delive[ l.t[ tourEval[i] ] - 1 ]){
						delive[ l.t[ tourEval[i] ] - 1] = true;
					 }
					 else{
						valido = false;
					 }
				  }
				  else{
					 valido = false;
				  }
				}
				for(int j = 1; j < l.k; j++){
				  if(!picked[j]) valido = false;
				  if(!delive[j]) valido = false;
				}

				// (3)
				double load = 0.0;
				for(int v = 0; v < (int)tourEval.size(); v++){
				  if( l.t[ tourEval[v] ] > 0 ){
					 load = load - l.items[ l.t[ tourEval[v] ] - 1 ].w;
				  }
				  if( l.s[ tourEval[v] ] > 0 ){
					 load = load + l.items[ l.s[ tourEval[v] ] - 1 ].w;
				  }
				  if(load < (-1)*MY_EPS) valido = false;
				  if(load > l.capacity) valido = false;
				}
				if(load > MY_EPS){
				  valido = false;
				}
				
				if(valido){
					double calcCost = 0.0;
					for(int i = 0; i < (int)tourEval.size(); i++){
						for(OutArcIt o(l.g, tourEval[i]); o != INVALID; ++o){
							if(l.g.target(o) == tourEval[(i+1) % (int)tourEval.size()]){
								calcCost += l.weight[o];
								break;
							}
						}
					}
					   
					if(calcCost < s.cost){
						s.tour.clear();
						unsigned int k;
						for(k = 0; k < tourEval.size(); k++){
							s.tour.push_back(tourEval[k]);
						}
						s.cost = calcCost;
					}
				}
			}
		}
		
		
		loopClock = clock();
		x++;
	}while ( ((double(loopClock - begin) / CLOCKS_PER_SEC) < tl) || (x < (int) l.items.size()) );
	
	s.lowerBound = s.cost;
	
	if((double(loopClock - begin) / CLOCKS_PER_SEC) < tl){
		s.upperBound = s.cost;
		return true;
	}
	
	return false;
}
//------------------------------------------------------------------------------
bool exact(const LpdTspInstance &l, LpdTspSolution  &s, int tl)
/* Implemente esta função, entretanto, não altere sua assinatura */
{
	
	return naive(l, s, tl);
}
//------------------------------------------------------------------------------
bool naive(const LpdTspInstance &instance, LpdTspSolution  &sol, int tl)
/*
 * Algoritmo ingênuo para o LPD-TSP. Ideia:
 * constrNaiveHeur(l, s)
 *    s.tour.push_back(l.depot)
 *    while(s.tour.size() < 2*l.k+1)
 *       v = argmin_{v' in V} {d_{(v,v')} | (v' é adj a v) e ((v' é s) ou (v' é t de i cujo s é u em l.tour))}
 *       l.tour.push_back(v)
 */
{
   DNode v,
         vl;

   double vval,
          vlval;

   int i;

   sol.tour.clear();
   sol.cost = 0.0;

   v = instance.depot;
   sol.tour.push_back(v);

   while((int)sol.tour.size() < 2 * instance.k + 1 && v != INVALID){
      v    = INVALID;
      vval = DBL_MAX;

      for(OutArcIt o(instance.g, sol.tour.back()); o != INVALID; ++o){
         vl    = instance.g.target(o);
         vlval = DBL_MAX;

         i = 0;
         while(i < (int)sol.tour.size() && vl != sol.tour[i]) i++;
         if(i < (int)sol.tour.size()) continue;

         if(instance.s[vl] > 0){
            vlval = instance.weight[o];
         }
         else if(instance.t[vl] > 0){
            i = 0;
            while(i < (int)sol.tour.size() && instance.t[vl] != instance.s[sol.tour[i]]){
               i++;
            }
            if(i < (int)sol.tour.size()){
               vlval = instance.weight[o];
            }
         }

         if(vlval < vval){
            v    = vl;
            vval = vlval;
         }
      }

      if(v != INVALID){
         sol.tour.push_back(v);
         sol.cost += vval;
      }
   }

   if(v == INVALID){
      sol.cost = DBL_MAX;
   }
   else{
      OutArcIt o(instance.g, sol.tour.back());
      for(; o != INVALID; ++o){
         if(instance.g.target(o) == sol.tour.front()) break;
      }
      if(o != INVALID){
         sol.cost += instance.weight[o];
      }
   }
	
	return false;
}
//------------------------------------------------------------------------------

