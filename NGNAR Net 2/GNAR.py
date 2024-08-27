import numpy as np

class NetworkModel:
    def __init__(self, alpha_order, beta_order, intercept=True, global_intercept=False):
        self.beta_order = beta_order
        self.intercept = intercept
        self.global_intercept = global_intercept
        if isinstance(alpha_order,int):
            self.alpha_order = alpha_order
            self.alpha_orders = [1 for _ in range(alpha_order)]
        else:
            self.alpha_order = len(alpha_order)
            self.alpha_orders = alpha_order # alpha_orders is [1,0,0,1,...]; alpha order is max time lag
        
    """
    def transformVTS(self, vts):
        # shape of X: n*p, col_1:x_t-1; col_2:2+beta_order_1: weighted t-1 stage-r neighbours, ...
        X = np.zeros(((len(vts)-self.alpha_order)*self.network.size,self.coef_order))
        for t in range(1,self.alpha_order+1):
            #alpha_t
            if self.coef_index[2*t-2] != 0 and self.coef_index[2*t-2] != self.coef_index[2*t-3]:
                X[:,self.coef_index[2*t-2]-1] = vts[self.alpha_order-t:-t,:].flatten("F")
            for i in range(self.coef_index[2*t-2], self.coef_index[2*t-1]):
                #each beta_t
                stage = i - self.coef_index[2*t-2]
                #for each node, compute weighted t-1 stage-r neighbours, and concat them to one vertical column
                X[:,i] = np.concatenate([vts[self.alpha_order-t:-t,self.network.neighbours[node][stage]].
                                         dot(self.network.neighbour_weights[node][stage])
                                       if len(self.network.neighbours[node]) > stage
                                       else [0 for i in range(vts.shape[0]-alpha_order)]
                                       for node in self.network.nodes])
        if self.intercept:
            if self.global_intercept:
                X = np.hstack((X,np.ones((X.shape[0],1))))
            else:
                X = np.hstack((X,np.repeat(np.identity(self.network.size), len(vts)-self.alpha_order, axis=0)))
        return X
        """

    def transformVTS(self,vts):
        X = []
        for t in range(self.alpha_order):
            X.append(np.transpose(vts[self.alpha_order-t-1:-t-1,:]@self.network.w_mats[1-self.alpha_orders[t]:1+self.beta_order[t]],axes=(1,2,0)).reshape(
                (-1,self.alpha_orders[t]+self.beta_order[t]),order="F")) #optimize this
        X = np.hstack(X)
        if self.intercept:
            if self.global_intercept:
                X = np.hstack((X,np.ones((X.shape[0],1))))
            else:
                X = np.hstack((X,np.repeat(np.identity(self.network.size), len(vts)-self.alpha_order, axis=0)))
        return X
    
    def return_coef(self, t, type="all"):
        if type=="alpha":
            if self.coef_index[2*t-2] != 0 and self.coef_index[2*t-2] != self.coef_index[2*t-3]:
                return self.coefs[self.coef_index[2*t-2]-1]
            else:
                return 0
        elif type == "beta":
            return self.coefs[self.coef_index[2*t-2]:self.coef_index[2*t-1]]
        elif type == "all":
            return self.coefs[self.coef_index[2*t-2]-1:self.coef_index[2*t-1]]

    def AIC(self):
        return 2*len(self.coefs)-2*self.loglik()
    def BIC(self):
        return len(self.coefs)*np.log(len(self.y))-2*self.loglik()
    def group_NTS(self, network, vts, size):
        self.network = network
        X = self.transformVTS(vts)
        T,n = len(vts),self.network.size
        self.k = int(np.ceil((T-self.alpha_order)/size))
        G = np.zeros(((T-self.alpha_order)*n,self.k*X.shape[1]))
        for i in range(self.k-1):
            for j in range(i+1):
                G[i*size*n:(i+1)*size*n,j*X.shape[1]:(j+1)*X.shape[1]] = X[i*size*n:(i+1)*size*n,:] 
        for j in range(self.k):
            G[(self.k-1)*size*n:,j*X.shape[1]:(j+1)*X.shape[1]] = X[(self.k-1)*size*n:,:]
        return G

from scipy.optimize import lsq_linear
from scipy.stats import multivariate_normal, norm
from itertools import chain, combinations


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def min_part(xs, r):
    xs = np.sort(xs)
    part = []
    ys = []
    for x in xs:
        for y in ys:
            if abs(x-y) >= r:
                part.append(ys)
                ys = []
                break
        ys.append(x)
    part.append(ys)
    return part
            

class GNAR(NetworkModel):
    def __init__(self, alpha_order, beta_order, intercept=True, global_intercept=False):
        super().__init__(alpha_order, beta_order,intercept,global_intercept)
        self.coef_index = np.cumsum([[self.alpha_orders[i]]+[beta_order[i]] for i in range(self.alpha_order)])
        self.coef_order = self.coef_index[-1]
        
    def simulate(self, network, initial_vts, length, coefs, error_cov_mat):
        self.network = network
        l = len(initial_vts)
        vts_sim = np.zeros((length+l,self.network.size))
        vts_sim[:l,:] = initial_vts
        for i in range(length):
            vts_sim[l+i,:] = self.transformVTS(vts_sim[l+i-self.alpha_order:l+i+1,:])@coefs + multivariate_normal(
                cov=error_cov_mat).rvs(1)
        return vts_sim[len(initial_vts):]
    
    def loglik(self):
        return np.sum(norm.logpdf(self.y,loc=self.X@self.coefs,scale=self.sigma2**.5))
    
    def r2(self):
        return r2_score(self.y, self.X@self.coefs)

    def adj_r2(self):
        return 1-((1-self.r2())*(len(self.y)-1)/(len(self.y)-len(self.coefs)-1))
        
    def fit(self, network, vts):
        # the order of cols of vts must match the order of network nodes
        # i.e. the first vts col is the ts for the first node in the network, etc
        self.network = network
        self.vts = vts
        self.vts_end = vts[-self.alpha_order:,:]
        self.y = vts[self.alpha_order:,:].flatten("F")
        self.X = self.transformVTS(vts)
        self.coefs,self.res = np.linalg.lstsq(self.X, self.y, rcond=None)[:2]
        self.vts_fitted = np.reshape((self.X@self.coefs),(-1,network.size),"F")
        #self.sigma2 = np.sum((self.y-self.X@self.coefs)**2)/(len(self.y)-len(self.coefs)) #assuming equal error var
        
    def compute_cov(self,X,y,coefs):
        T = int(len(y)/self.network.size)
        residual = np.reshape((y-X@coefs),(-1,self.network.size),"F")
        res_cov_mat = residual.T@residual/T
        var_X_res = 0
        X = X[np.concatenate([np.array([i*T for i in range(self.network.size)])+t for t in range(T)])]
        y = y[np.concatenate([np.array([i*T for i in range(self.network.size)])+t for t in range(T)])]
        for i in range(self.network.size):
            for j in range(self.network.size):
                xti_xtj = np.mean([np.outer(X[i+k*self.network.size],X[j+k*self.network.size]) for k in range(T)],axis=0)
                var_X_res += res_cov_mat[i,j]*xti_xtj
        X_cov_mat = X.T@X/T
        inv_X_cov_mat = np.linalg.inv(X.T@X)*T
        coef_cov_mat = inv_X_cov_mat@var_X_res@inv_X_cov_mat/T
        inv_coef_cov_mat = X_cov_mat@np.linalg.solve(var_X_res,X_cov_mat)*T
        return coef_cov_mat,inv_coef_cov_mat
        
    def predict(self, length, nodes=None, vts_end =None):
        if vts_end is None:
            vts_end = self.vts_end
        vts_pred = np.zeros((length+self.alpha_order,self.network.size))
        vts_pred[:self.alpha_order,:] = vts_end
        for i in range(length):
            vts_pred[self.alpha_order+i,:] = self.transformVTS(vts_pred[i:self.alpha_order+i+1,:])@self.coefs
        if nodes is None:
            return vts_pred[self.alpha_order:,:]
        else:
            return vts_pred[self.alpha_order:,nodes]


    
    def return_coef(self, t, type="all"):
        if type=="alpha":
            if self.coef_index[2*t-2] != 0 and self.coef_index[2*t-2] != self.coef_index[2*t-3]:
                return self.coefs[self.coef_index[2*t-2]-1]
            else:
                return 0
        elif type == "beta":
            return self.coefs[self.coef_index[2*t-2]:self.coef_index[2*t-1]]
        elif type == "all":
            return self.coefs[self.coef_index[2*t-2]-1:self.coef_index[2*t-1]]
        
    def initial_fit(self,network,vts,size,level=0.05):
        self.network = network
        self.vts = vts
        self.size = size
        self.T = len(self.vts)-self.alpha_order
        self.y = vts[self.alpha_order:,:].flatten("F")
        self.X = self.transformVTS(vts)
        self.block_number = int(np.ceil(self.T/size))
        index = np.concatenate([np.arange(size)+i*self.T for i in range(self.network.size)])
        coefs_prev = lsq_linear(self.X[index], self.y[index]).x
        self.cpts_candidate = []
        for l in range(1,self.block_number):
            index = np.concatenate([np.arange(l*size,min(self.T,(l+1)*size))+i*self.T for i in range(self.network.size)])
            coefs_new = lsq_linear(self.X[index], self.y[index]).x
            coef_cov_mat,inv_coef_cov_mat = self.compute_cov(self.X[index], self.y[index],coefs_new)
            test_chi2 = np.dot((coefs_new-coefs_prev),(inv_coef_cov_mat@(coefs_new-coefs_prev)))
            p = chi2(df=len(coefs_prev)).sf(test_chi2)
            #print(p)
            #print(coef_cov_mat)
            #print(2*norm().sf(abs(coefs_new-coefs_prev)/np.array([coef_cov_mat[i,i] for i in range(len(coefs_new))])**(1/2)))
            if p<level:
                coefs_prev = coefs_new
                self.cpts_candidate.append(l*size)
        #print(self.cpts_candidate)
            
    def minimize(self,X,y):
        return lsq_linear(X,y).cost
            
    def LIC(self,cpts_subset,window,omega):
        #local information criterion
        rss = 0
        n = 0
        for cpt in self.cpts_candidate:
            if (cpt in cpts_subset):
                index_1 = np.concatenate([np.arange(cpt-window,cpt)+i*(self.T) for i in range(self.network.size)])
                index_2 = np.concatenate([np.arange(cpt,min(self.T,cpt+window))+i*(self.T) for i in range(self.network.size)])
                rss += (lsq_linear(self.X[index_1],self.y[index_1]).cost+
                lsq_linear(self.X[index_2],self.y[index_2]).cost)
                n += len(index_1) + len(index_2)
            else:
                index = np.concatenate([np.arange(cpt-window,min(self.T,cpt+window))+i*(self.T) for i in range(self.network.size)])
                rss += lsq_linear(self.X[index],self.y[index]).cost
                n += len(index)
        if omega == "AIC":
            lic = n*np.log(rss/n) + len(cpts_subset)*self.X.shape[1]*2
        elif omega == "BIC":
            lic = n*np.log(rss/n) + len(cpts_subset)*self.X.shape[1]*np.log(n)
        else:
            lic = n*np.log(rss/n) + len(cpts_subset)*self.X.shape[1]*omega
        return lic
            
    def local_screening(self,window,omega):
        lowest_lic = np.inf
        for cpts_subset in powerset(self.cpts_candidate):
            lic = self.LIC(cpts_subset,window,omega)
            if lic < lowest_lic:
                cpts = cpts_subset
                lowest_lic = lic
        self.cpts = cpts
    
    
    def exhaustive_search(self,window):
        partition = min_part(self.cpts,2*window)
        cpts = []
        for ts in partition:
            if len(ts) == 1:
                l,u = ts[0]-window, min(self.T,ts[0]+window)
            else:
                l,u = min(ts),max(ts)
            lowest_cost = np.inf
            for t in range(l,u):
                index_1 = np.concatenate([np.arange(l,t)+i*(len(self.vts)-self.alpha_order) for i in range(self.network.size)])
                index_2 = np.concatenate([np.arange(t,u)+i*(len(self.vts)-self.alpha_order) for i in range(self.network.size)])
                cost = lsq_linear(self.X[index_1],self.y[index_1]).cost + lsq_linear(self.X[index_2],self.y[index_2]).cost
                if cost < lowest_cost:
                    best_t = t
                    lowest_cost = cost
            cpts.append(best_t)
        self.cpts = np.array(cpts)+self.alpha_order
        
    def cpts_detect(self,network,vts,size,window,omega):
        self.initial_fit(network,vts,size)
        if len(self.cpts_candidate) == 0:
            return []
        self.local_screening(window, omega)
        if len(self.cpts) == 0:
            return []
        self.exhaustive_search(window)
        return self.cpts
            

def GNAR_sim_piecewise(networks, alpha_order, beta_order, coefs_list, intercept,global_intercept,error_cov_mat_list, length_list, burn_in = 100):
    num_piece = len(coefs_list)
    vts = np.zeros((sum(length_list)+burn_in+alpha_order,networks[0].size))
    vts[:alpha_order] = norm(loc=0).rvs(size=(alpha_order,networks[0].size))
    current_index = alpha_order
    length_list[0] += burn_in
    for l in range(num_piece):
        model = GNAR(alpha_order,beta_order,intercept, global_intercept)
        vts[current_index:current_index+length_list[l]]=model.simulate(networks[l],vts[:current_index],length_list[l],coefs_list[l],error_cov_mat_list[l])
        current_index = current_index+length_list[l]
    return vts[burn_in+alpha_order:]


def softplus(x,c=1,limit=30):
    sp = np.zeros(x.shape)
    sp[x/c>limit] = x[x>limit]
    sp[x/c<=limit] = c*np.log1p(np.exp(x[x<=limit]/c))
    return sp

import tensorflow as tf
import tensorflow_probability as tfp

def norm_tf(x, order = 1):
    return tf.math.pow(tf.reduce_sum(tf.math.pow(tf.abs(x),order)),1/order)


class NGNAR_TF(NetworkModel):
    def __init__(self, alpha_order, beta_order, intercept=True, global_intercept=True,link="softplus",seasonal=False):
        self.beta_order = beta_order
        self.intercept = intercept
        self.global_intercept = global_intercept
        
        if isinstance(alpha_order,int):
            self.alpha_order = alpha_order
            self.alpha_orders = [1 for _ in range(alpha_order)]
        else:
            self.alpha_order = len(alpha_order)
            self.alpha_orders = alpha_order # alpha_orders is [1,0,0,1,...]; alpha order is max time lag
        
        self.coef_index = np.cumsum([[self.alpha_orders[i]]+[beta_order[i]] for i in range(self.alpha_order)])
        self.coef_order = self.coef_index[-1]
            
        self.seasonal = seasonal #seasonal adjustment 
        #self.inv =inv # controls whether the link function is applied to the vts when fitting
        if link == "softplus":
            self.link = lambda x:tf.math.softplus(x)
            self.link_np = lambda x: softplus(x)
        elif link == "identity":
            self.link = lambda x:x
            self.link_np = self.link
        elif link == "relu":
            self.link = tf.nn.relu
            self.link_np = relu
    
    def get_fitted(self):
        y = self.link(tf.linalg.matmul(self.X,self.coefs))
        if self.seasonal:
            return y + self.season_seq(self.t)
        else:
            return y
        
    def season_seq(self,t):
        # the last term is to keep the result non-neg
        return self.coefs_season[0]*tf.sin(self.coefs_season[1]*t-self.coefs_season[2])+self.coefs_season[0] 
    
    def MSE(self):
        return tf.keras.losses.MeanSquaredError()(self.y, self.get_fitted())+self.coef_penalty()
    
    def r2(self):
        return r2_score(self.y, self.get_fitted())

    def adj_r2(self):
        return 1-((1-self.r2())*(len(self.y)-1)/(len(self.y)-len(self.coefs)-1))
    
    def loglik(self):
        return np.sum(poisson.logpmf(self.y,mu=self.y_fitted))

    
    def neg_loglik(self):
        return -tf.math.reduce_sum(tfp.distributions.Poisson(rate=self.get_fitted()).log_prob(self.y))+self.coef_penalty()
    
    def fit_MLE(self, network, vts):
        """
        fit the NGNAR model using IWLS, 
        only poisson distribution is implemented 
        """
        glm = tfp.glm.CustomExponentialFamily(tfp.distributions.Poisson, tf.math.softplus)
        self.network = network
        
        self.vts_end = vts[-self.alpha_order:,:]
        self.y = tf.convert_to_tensor(vts[self.alpha_order:,:].flatten("F"),dtype=np.float64)
        self.X = tf.convert_to_tensor(self.transformVTS(vts),dtype=np.float64)
        self.coefs = tf.Variable(np.linalg.lstsq(self.X, self.y,None)[0])
        self.coefs = tfp.glm.fit(self.X,self.y,glm,self.coefs,fast_unsafe_numerics=False)[0].numpy()
    
    def fit(self, network, vts, loss = "MSE", reg_order=None, reg_coef=None):
        """
        fit the NGNAR model by minimizing a loss function,
        loss=neg_loklik => mle;
        loss=mse => ls
        """
        # the number of cols of vts must match the size of the network
        # i.e. the first vts col is the ts for the first node in the network, etc
        self.network = network
        self.vts_end = vts[-self.alpha_order:,:]
        self.y = tf.convert_to_tensor(vts[self.alpha_order:,:].flatten("F").reshape(-1,1))
        self.X = tf.convert_to_tensor(self.transformVTS(vts))
        self.coefs = tf.Variable(np.linalg.lstsq(self.X, self.y,None)[0].reshape(-1,1))
        
        self.t = tf.tile(tf.reshape(tf.range(self.alpha_order,len(vts),dtype=tf.float64),(-1,1)),
                         tf.constant([network.size,1]))
        self.coefs_season = tf.Variable(np.array([1.,1.,0.]))#0-scale, 1-time 2-transform
        
        if loss == "MSE":
            loss = self.MSE
        elif loss == "neg_loglik":
            loss = self.neg_loglik
            
        self.coef_penalty=lambda:0
        
        if reg_order is not None:
            self.coef_penalty = lambda : self.reg_coef* tf.norm(self.coefs[:self.coef_order],reg_order)
            
            if reg_coef is None:
                self.reg_coef = 1
            else:
                self.reg_coef = reg_coef
            
        #trace_fn = lambda traceable_quantities: {'x': self.coefs}
        trace = tfp.math.minimize(loss, num_steps=10000,
                                  optimizer=tf.optimizers.Adam(),
                                  convergence_criterion=(tfp.optimizer.convergence_criteria.LossNotDecreasing(atol=0.01)))
        self.coefs = self.coefs.numpy()
        self.y_fitted = self.get_fitted().numpy()
        #self.y_fitted[self.y_fitted==0]+=0.001
        self.vts_fitted = np.reshape(self.link(self.X@self.coefs),(-1,network.size),"F")
        
    def predict(self, length, nodes=None):
        vts_pred = np.zeros((length+self.alpha_order,self.network.size))
        vts_pred[:self.alpha_order,:] = self.vts_end
        for i in range(length):
            vts_pred[self.alpha_order+i,:] = self.link_np(
                np.log(self.transformVTS(vts_pred[i:self.alpha_order+i+1,:]+1))@self.coefs).reshape(-1)
        if self.seasonal:
            """
            wrong!! relations should be related to X_t itself!
            """
            return vts_pred[self.alpha_order:,:]+self.season_seq(np.linspace(len(self.t),len(self.t)+length,1))
        else:
            return vts_pred[self.alpha_order:,:]
    def return_coef(self, t, type="all"):
        if type=="alpha":
            if self.coef_index[2*t-2] != 0 and self.coef_index[2*t-2] != self.coef_index[2*t-3]:
                return self.coefs[self.coef_index[2*t-2]-1]
            else:
                return 0
        elif type == "beta":
            return self.coefs[self.coef_index[2*t-2]:self.coef_index[2*t-1]]
        elif type == "all":
            return self.coefs[self.coef_index[2*t-2]-1:self.coef_index[2*t-1]]

def simulate_NGNAR(network, alpha_order, beta_order, intercept, global_intercept, length, coefs, link="softplus", burn_in=100):
    vts = np.zeros((length+burn_in+alpha_order,network.size))
    vts[:alpha_order] = np.random.poisson(size=(alpha_order,network.size))
    model = NGNAR_TF(alpha_order,beta_order,intercept,global_intercept,link)
    model.network = network
    for i in range(alpha_order, length+burn_in+alpha_order):
        vts[i] = np.random.poisson( model.link_np(model.transformVTS(vts[i-alpha_order:i+1,:])@coefs) ).reshape(-1)
    return vts[burn_in+alpha_order:]