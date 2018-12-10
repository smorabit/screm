import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import time
import matplotlib
from multiprocessing import Pool
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from matplotlib import rc

class Screm:

	#empty constructor
	def __init__(self):
		pass

	def delta_t(self, l):
	    return -math.log(np.random.uniform())/l

	def simulate_reactions(self, v, reaction_net, iterations, verbose=False):
    
	    start_time = time.time()
	    np.random.seed()
	    A = reaction_net["stoichiometry"]
	    reactants = reaction_net["reactants"]
	    reaction_rates = reaction_net["rates"]
	    reaction_relations = reaction_net["relations"]
	    
	    times = [0]
	    
	    for i in range(iterations):
	        
	        #get relevant reactants
	        reactant_indices = [j for j, item in enumerate(np.ravel(v[-1])) if item !=0]
	        cur_reactants = [reactants[j] for j in reactant_indices]
	        cur_reactant_amounts = {r: v[-1,j] for j, r in enumerate(reactants)}
	        
	        #compute possible next reactions
	        next_reactions = [rxn for rxn, req in reaction_relations.items() if set(req).issubset(set(cur_reactants))]
	        
	        #exit if we have there are no possible next reactions:
	        if len(next_reactions) == 0:
	            break
	                
	        #compute rxn coefficients for possible next reactions:
	        rxn_rates = {rxn: reaction_rates[rxn]*math.exp(sum(map(math.log,[v[-1,reactants.index(r)] for r in reaction_relations[rxn]]))) for rxn in next_reactions}
	        
	        #compute delta ts:
	        dt = {rxn: self.delta_t(rate) for rxn, rate in rxn_rates.items()}
	        
	        #find fastest dt:
	        while 1:
	            next_state = min(dt, key=dt.get)
	            update_v = v[-1]+A[:,next_state].T
	            if -1 not in np.ravel(update_v):
	                break
	            else:
	                del dt[next_state]
	        
	        #update times:
	        times.append(times[-1] + dt[next_state])
	        
	        #update v:
	        v = np.vstack((v, update_v))
	        
	        if verbose == True:
	            print("v[-1]:", np.ravel(v[-1]))
	            print("cur reactants", cur_reactants)
	            print("cur_reactant_amounts", cur_reactant_amounts)
	            print("next reactions", next_reactions)
	            print("rxn_rates", rxn_rates)
	            print("dt:", dt)
	            print("fastest dt:", min(dt, key=dt.get))
	            print()
	    end_time = time.time()
	    #print("time elapsed: {}".format(end_time-start_time))
	    return times, v, i+1

	def plot_simulation(self, T,V,I, reaction_net, metabolites="all", figsize=(7,3), dpi=170):
	    
	    A = reaction_net["stoichiometry"]
	    reactants = reaction_net["reactants"]
	    reaction_rates = reaction_net["rates"]
	    reaction_relations = reaction_net["relations"]
	    
	    fig, axes = plt.subplots(1,1,figsize=figsize, dpi=dpi)

	    cm = plt.get_cmap('jet')
	    color_dict = {r: cm(1.*j/A.shape[0]) for j, r in enumerate(reactants)}
	    
	    #which reactants do we want to plot?
	    if metabolites != "all":
	        plot_metabolites = [reactants.index(m) for m in metabolites]
	    else:
	        plot_metabolites = [i for i in range(A.shape[0])]
	        metabolites = [r for r in reactants]
	    
	    for j in plot_metabolites:
	        axes.plot(T, np.ravel(V[:,j]), color=color_dict[reactants[j]])
	    axes.set_title("{} iterations".format(I))
	    axes.spines['top'].set_visible(False)
	    axes.spines['right'].set_visible(False)
	    axes.set_xticks([])
	    axes.set_xlabel("time")
	    axes.set_ylabel("molecules")
	    axes.legend(metabolites, loc='center left', bbox_to_anchor=(1, 0.5))
	    plt.tight_layout()
	    
	def plot_simulations_violin(self, model_results, init, ax, datatype="iter", poly=1):
	    
	    if datatype == "iter":
	        data = [[I for T,V,I in model_results[i]] for i in range(len(init))]
	        means = [np.mean(i) for i in data]
	    elif datatype == "time":
	        data = [[T[-1] for T,V,I in model_results[i]] for i in range(len(init))]
	        means = [np.mean(t) for t in data]

	    violin = ax.violinplot(data, vert=True, showmeans=True);
	    for i, v in enumerate(violin['bodies']):
	        v.set_edgecolor('k')
	        v.set_alpha(0.7)

	    # change the line color from blue to black
	    for partname in ('cbars','cmins','cmaxes', 'cmeans'):
	        vp = violin[partname]
	        vp.set_edgecolor('k')
	        vp.set_linewidth(1)

	    #draw best fit linear line:
	    z, residuals, _, _, _ = np.polyfit([i for i in range(1, len(init)+1)], means, poly, full=True)
	    fit = np.polyval(z, [i for i in range(1,len(init)+1)])
	    ax.plot([i for i in range(1, len(init_glucose)+1)], fit)
	    r2 = r2_score(fit, means)

	    #format spines:
	    ax.spines['top'].set_visible(False)
	    ax.spines['right'].set_visible(False)

	    #label things:
	    if datatype == "iter":
	        ax.set_ylabel("iterations to convergence")
	    elif datatype == "time":
	        ax.set_ylabel("time (in silico)")
	    ax.set_yticks([])
	    ax.set_xlabel("initial glucose molecules")
	    ax.set_xticks([2*i-1 for i in range(1,int(len(init)/2)+1)])
	    ax.set_xticklabels([g for i, g in enumerate(init_glucose) if i % 2 == 0])
	    ax.annotate("r-squared:{}".format(str(r2)[:5]),
	                xy=(1, 0), xycoords='axes fraction',
	                xytext=(-20, 20), textcoords='offset pixels',
	                horizontalalignment='right',
	                verticalalignment='bottom')

	def plot_residuals(self, model_results, init, ax, datatype="iter", poly=1):
	    
	    if datatype == "iter":
	        data = [[I for T,V,I in model_results[i]] for i in range(len(init))]
	        means = [np.mean(i) for i in data]
	    elif datatype == "time":
	        data = [[T[-1] for T,V,I in model_results[i]] for i in range(len(init))]
	        means = [np.mean(t) for t in data]
	    
	    #compute best fit line:
	    z, residuals, _, _, _ = np.polyfit([i for i in range(1, len(init)+1)], means, poly, full=True)
	    fit = np.polyval(z, [i for i in range(1,len(init)+1)])
	    r2 = r2_score(fit, means)
	    
	    ax.scatter([i for i in range(1,len(init)+1)],fit - means)

	    #format:
	    ax.spines['top'].set_visible(False)
	    ax.spines['right'].set_visible(False)
	    ax.set_ylabel("residuals")
	    ax.set_yticks([])
	    ax.set_xlabel("initial glucose molecules")
	    ax.set_xticks([2*i-1 for i in range(1,int(len(init)/2)+1)])
	    ax.set_xticklabels([g for i, g in enumerate(init) if i % 2 == 0])
	   
	def reactant_start_end_times(self, model_results, reaction_net):
	    start_end_list = []
	    for run in model_results:
	        start_end = {r:[[-1 for x in run],[-1 for x in run]] for r in reaction_net["reactants"]}
	        for k, stuff in enumerate(run):
	            T,V,I = stuff
	            for i in range(I):
	                for j, r in enumerate(reaction_net["reactants"]):

	                    #did this reactant appear for the first time?
	                    if V[i,j] > 0 and start_end[r][0][k] == -1:
	                        start_end[r][0][k] = T[i]

	                    #did this reactant disappear?
	                    elif V[i,j] == 0 and start_end[r][0][k] > -1 and start_end[r][1][k] == -1:
	                        start_end[r][1][k] = T[i]

	        start_end_list.append(start_end)
	    return start_end_list

	def plot_reactant_start_end_times(self, sel, reactant, ax, end=True):
	    
	    if end == True:
	        data = [se[reactant][1] for se in sel]
	    else:
	        data = [se[reactant][0] for se in sel]

	    violin = ax.violinplot(data, vert=True, showmeans=True);
	    for i, v in enumerate(violin['bodies']):
	        v.set_edgecolor('k')
	        v.set_alpha(0.7)

	    # change the line color from blue to black
	    for partname in ('cbars','cmins','cmaxes', 'cmeans'):
	        vp = violin[partname]
	        vp.set_edgecolor('k')
	        vp.set_linewidth(1)
	    
	    ax.spines['top'].set_visible(False)
	    ax.spines['right'].set_visible(False)
