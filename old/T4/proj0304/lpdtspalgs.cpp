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
#include <gurobi_c++.h>




bool naive(const LpdTspInstance &l, LpdTspSolution  &s, int tl);




//------------------------------------------------------------------------------
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
	clock_t begin = clock();
	
	DNode deposito = l.depot;
	double capacity = l.capacity;
	
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
	
	// criando variaveis
	vector<GRBVar> x(l.m); // indica se a aresta pertence a solucao
	vector<GRBVar> f(l.m); // indica o peso carregado ao passar na aresta
	GRBVar total;
	
	for(ListDigraph::ArcIt e(l.g); e!=INVALID; ++e){
		x[l.g.id(e)] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		f[l.g.id(e)] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	}
	total = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
	
	model.update();
	
	// criando objetivo
	GRBLinExpr expr1;
	for(ListDigraph::ArcIt e(l.g); e!=INVALID; ++e){
		expr1 += (l.weight[e] * x[l.g.id(e)]);
	}
	model.setObjective(expr1, GRB_MINIMIZE);
	
	model.update();
	
	// criando restricoes
	model.addConstr(expr1 == total); // somatorio(peso[a] * x[a]) == total
	
	for(ListDigraph::NodeIt v(l.g); v!=INVALID; ++v){
		GRBLinExpr expr2;
		GRBLinExpr expr3;
		GRBLinExpr expr4;
		for(OutArcIt e(l.g, v); e != INVALID; ++e){
			expr2 += x[l.g.id(e)];
			expr4 += f[l.g.id(e)];
		}
		for(InArcIt e(l.g, v); e != INVALID; ++e){
			expr3 += x[l.g.id(e)];
			expr4 -= f[l.g.id(e)];
		}
		model.addConstr(expr2 == 1.0); // somatorio(x[a]) das arestas que saem == 1
		model.addConstr(expr3 == 1.0); // somatorio(x[a]) das arestas incidentes == 1
		
		if(v!= deposito){
			if(l.s[v] > 0){
				// somatorio(f[a]) das arestas q saem - somatorio(f[b]) das arestas incidentes - item coletado = 0
				model.addConstr(expr4 - l.items[l.s[v] - 1].w == 0.0);
			}
			
			if(l.t[v] > 0){
				// somatorio(f[a]) das arestas q saem - somatorio(f[b]) das arestas incidentes + item entregue = 0
				model.addConstr(expr4 + l.items[l.t[v] - 1].w == 0.0);
			}
		}else{
			model.addConstr(expr4 == 0.0);
		}
	}
	
	for(ListDigraph::ArcIt e(l.g); e!=INVALID; ++e){
		GRBLinExpr expr5 = f[l.g.id(e)];
		model.addConstr(expr5 <= capacity); // f[a] <= capacidade do caixeiro
		model.addConstr(expr5 >= 0.0); // f[a] >= 0
		
		if(l.g.source(e) != deposito && l.g.target(e) != deposito){
			GRBLinExpr expr6 = f[l.g.id(e)];
			model.addConstr(expr6 >= x[l.g.id(e)]); // f[a] >= x[a]
		}else{ // arestas que ligam o deposito nao devem ter carga
			GRBLinExpr expr7 = f[l.g.id(e)];
			model.addConstr(expr7 == 0.0); // f[a] = 0
		}
	}
	
	model.update();
	
	model.getEnv().set(GRB_DoubleParam_TimeLimit, tl);
	
	bool feasible = false;
	GRBConstr lastConstr;
	GRBLinExpr lastExpr;
	int lastValue = -1;
	do{ // o resolvedor sera chamado quantas vezes forem necessarias para obter uma solucao viavel
		try {
			s.tour.clear();
			model.optimize();
			
			DNode atual = deposito;
			s.cost = total.get(GRB_DoubleAttr_X);
			s.lowerBound = s.cost;
			s.upperBound = s.cost;
			do{// atualiza o vetor s.tour
				s.tour.push_back(atual);
				
				for(OutArcIt e(l.g, atual); e != INVALID; ++e){
					if(x[l.g.id(e)].get(GRB_DoubleAttr_X) == 1.0){
						atual = l.g.target(e);
						break;
					}
				}
			}while(atual != deposito);

			feasible = true;

			for(int i = 0; i < l.k && feasible; i++){
				Item item = l.items[i];
				int fromIndex = -1;
				int toIndex = -1;
				
				// verifica a coerencia na entrega dos itens
				for(unsigned int j = 0; j < s.tour.size(); j++){
					ListDigraphBase::Node node = s.tour[j];
					if(node == item.s){
						fromIndex = j;
					}
					if(node == item.t){
						toIndex = j;
					}
				}
				
				if(fromIndex > toIndex && fromIndex != -1 && toIndex != -1){ // entregue antes de ser coletado
					feasible = false;
					
					for(ListDigraph::ArcIt e(l.g); e!=INVALID; ++e){
						if(l.g.target(e) == item.s && x[l.g.id(e)].get(GRB_DoubleAttr_X) == 1.0){
							GRBLinExpr expr8;
							expr8 += x[l.g.id(e)];
							lastConstr = model.addConstr(expr8 == 0.0);
							lastExpr = expr8;
							lastValue = 0;
							model.update();
						}
					}
				}else{
					if(fromIndex == -1){ // nunca coletado
						feasible = false;
						Arc e = INVALID;
						int i = 1;
						while(e == INVALID){
							e = findArc(l.g, s.tour[s.tour.size() - i], item.s);
							i++;
						}
						GRBLinExpr expr9 = x[l.g.id(e)];
						lastConstr = model.addConstr(expr9 == 1.0);
						lastExpr = expr9;
						lastValue = 1;
						model.update();
					}else{
						if(toIndex == -1){ // nunca entregue
							feasible = false;
							Arc e = INVALID;
							int i = 1;
							while(e == INVALID){
								e = findArc(l.g, s.tour[s.tour.size() - i], item.t);
								i++;
							}
							GRBLinExpr expr10 = x[l.g.id(e)];
							lastConstr = model.addConstr(expr10 == 1.0);
							lastExpr = expr10;
							lastValue = 1;
							model.update();
						}
					}
				}
			}
			
			clock_t posAnalise = clock();
	
			if( (double(posAnalise - begin) / CLOCKS_PER_SEC) > tl ){
				s.lowerBound = model.get(GRB_DoubleAttr_MinBound);
				s.upperBound = s.cost;
				return false;
			}
			
			// atualiza o tempo maximo de otimizacao
			model.getEnv().set(GRB_DoubleParam_TimeLimit, tl - ((posAnalise - begin) / CLOCKS_PER_SEC));
			
		}catch(GRBException e) {
			clock_t posAnalise = clock();
	
			if( (double(posAnalise - begin) / CLOCKS_PER_SEC) > tl ){
				s.lowerBound = model.get(GRB_DoubleAttr_MinBound);
				s.upperBound = s.cost;
				return false;
			}
			
			// atualiza o tempo maximo de otimizacao
			model.getEnv().set(GRB_DoubleParam_TimeLimit, tl - ((posAnalise - begin) / CLOCKS_PER_SEC));
			
			// retira a ultima restricao adicionada
			model.remove(lastConstr);
			
			if(lastValue == 0){
				model.addConstr(lastExpr == 1.0);
				lastValue = 1;
			}else{
				model.addConstr(lastExpr == 0.0);
				lastValue = 0;
			}
			
			model.update();
		}
	}while(!feasible);
	
	s.lowerBound = s.cost;
	s.upperBound = s.cost;
		
	return true;
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

         if(instance.s[vl] > 0){  // If DNode vl is start of an item
            vlval = instance.weight[o];
         }
         else if(instance.t[vl] > 0){  // If DNode vl is término of an item
            i = 0;
            while(i < (int)sol.tour.size() && instance.t[ vl ] != instance.s[ sol.tour[i] ]){  // Look for the start DNode of the item which terminates in DNode vl
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
