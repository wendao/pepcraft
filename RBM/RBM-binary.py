# -*- coding: utf-8 -*-

import sys, os, argparse
from math import sqrt, log
import mxnet as mx
import mxnet.ndarray as nd
import numpy as np

def int2xy(conf):
    L = len(conf)
    x = np.zeros([L+2, 2])
    x[0,0] = 0
    x[0,1] = 0
    x[1,0] = 1
    x[1,1] = 0
    for i in range(L):
        if conf[i] == 1:
            #forward
            x[i+2,0] = 2*x[i+1,0]-x[i,0]
            x[i+2,1] = 2*x[i+1,1]-x[i,1]
        elif conf[i] == 3:
            #back
            x[i+2,0] = x[i,0]
            x[i+2,1] = x[i,1]
        else:
            dx1 = x[i+1,0] - x[i,0]
            dy1 = x[i+1,1] - x[i,1]
            if conf[i] == 2:
                #right
                if dx1 == 0:
                    dy = 0
                    if dy1>0: dx = -1
                    else: dx = 1
                elif dy1 == 0:
                    dx = 0
                    if dx1>0: dy = 1
                    else: dy = -1
                else:
                    print("error!")
            elif conf[i] == 0:
                #left
                if dx1 == 0:
                    dy = 0
                    if dy1>0: dx = 1
                    else: dx = -1
                elif dy1 == 0:
                    dx = 0
                    if dx1>0: dy = -1
                    else: dy = 1
                else:
                    print("error!")
            x[i+2,0] = x[i+1,0] + dx
            x[i+2,1] = x[i+1,1] + dy
    return x

map_i2b = {'0':'00', '1':'01', '2':'10', '3':'11'}
def int2binary(int_conf):
    bin_conf = ""
    for i in range(len(int_conf)):
        bin_conf += map_i2b[int_conf[i]]
    return bin_conf

map_b2i = {'00':'0', '01':'1', '10':'2', '11':'3'}
def binary2int(bin_conf):
    int_conf = ""
    for i in range(len(bin_conf)/2):
        j = i * 2
        int_conf += map_b2i[bin_conf[j:j+2]]
    return int_conf

def cal_rg2(xy):
    cm = np.sum(xy, axis=0) / xy.shape[0]
    return np.mean(np.sum((xy-cm)**2, axis=1))

class RBM(object):
    def __init__(self, n_vis=2, n_hid=1, n_val=2, ctx=mx.cpu()):
        self.n_vis = n_vis  # num of units in visible (input) layer
        self.n_hid = n_hid  # num of units in hidden layer
        self.n_val = n_val  # num of values for each node
        self.n_node = n_vis * 2 # each visnode has n_val nodes

        self.ctx = ctx

        a = sqrt(6. / ( self.n_vis * n_val + n_hid ))

        self.W = nd.random_uniform(low=-a, high=a, shape=(self.n_node, n_hid), ctx=self.ctx)
        self.hb = nd.zeros([n_hid], ctx=self.ctx)        # initialize h bias 0
        self.vb = nd.zeros([self.n_node], ctx=self.ctx)  # initialize v bias 0
        self.dW = nd.zeros([self.n_node, n_hid], ctx=self.ctx)
        self.dh = nd.zeros([n_hid], ctx=self.ctx)
        self.dv = nd.zeros([self.n_node], ctx=self.ctx)

        #for KL
        self.enum_states = None
        self.prob_states = None
        self.prob_RGs = None

        self.L_coeff = 0.0
        self.M_coeff = 0.0

    def contrastive_divergence(self, input, lr=0.1, cdk=1, batch_size=None, shuffle=False):
        n_sample = input.shape[0]
        if batch_size == 0: batch_size = n_sample

        labels = nd.ones([n_sample, 1], ctx=self.ctx)
        dataiter = mx.io.NDArrayIter(input, labels, batch_size, shuffle, last_batch_handle='discard')

        for batch in dataiter:
            sub = batch.data[0]

            ph_prob, ph_sample = self.sample_h_given_v(sub)
            chain_start = ph_sample

            for step in range(cdk):
                if step == 0:
                    nv_prob, nv_sample, nh_prob, nh_sample = self.gibbs_hvh(chain_start)
                else:
                    nv_prob, nv_sample, nh_prob, nh_sample = self.gibbs_hvh(nh_sample)

            if self.M_coeff > 0:
                self.dW *= self.M_coeff
                self.dv *= self.M_coeff
                self.dh *= self.M_coeff
                self.dW += (nd.dot(sub.T, ph_prob) - nd.dot(nv_sample.T, nh_prob)) * lr / batch_size
                self.dv += nd.mean(sub - nv_sample, axis=0) * lr
                self.dh += nd.mean(ph_prob - nh_prob, axis=0) * lr
            else:
                self.dW = (nd.dot(sub.T, ph_prob) - nd.dot(nv_sample.T, nh_prob)) * lr / batch_size
                self.dv = nd.mean(sub - nv_sample, axis=0) * lr
                self.dh = nd.mean(ph_prob - nh_prob, axis=0) * lr

            self.W  += self.dW
            self.vb += self.dv
            self.hb += self.dh

            self.W_decay(lr)
        return

    def W_decay(self, lr):
        #go through weights
        if self.L_coeff>0: #L2
            self.W -= self.L_coeff * lr * self.W
        elif self.L_coeff<0: #L1
            self.W += self.L_coeff * lr * nd.sign(self.W)
        else:
            #fix boundary
            self.W = nd.clip(self.W, -10.0, 10.0)
        return

    def check_status(self, input, epoch):
        n_sample = input.shape[0]

        ph_prob, ph_sample = self.sample_h_given_v(input)
        nv_prob, nv_sample, nh_prob, nh_sample = self.gibbs_hvh(ph_sample)
        error = nd.sum((input - nv_sample) ** 2) / n_sample
        #use logsoftmax if nan
        cross = -nd.mean(nd.sum(input * nd.log(nv_prob+1e-10), axis=1))
        freeE = self.get_free_energy(input)

        sys.stdout.write( "Training: " )
        sys.stdout.write( "epoch= %d " % epoch )
        sys.stdout.write( "cross= %f " % cross.asnumpy()[0] )
        sys.stdout.write( "error= %f " % error.asnumpy()[0] )
        sys.stdout.write( "freeE= %f " % freeE.asnumpy()[0] )

        if self.enum_states is not None:
            sys.stdout.write( "KL= %f " % self.check_KL() )
        if self.prob_RGs is not None:
            sys.stdout.write( "rgKL= %f " % self.check_rgKL(nv_sample) )

        sys.stdout.write("\n")
        return

    def load_enum_states(self, fn):
        lines = open(fn, 'r').readlines()
        n_states = int(lines[0])
        dat_lst = []

        self.prob_states = nd.zeros([n_states], ctx=self.ctx)

        for i in range(1, n_states+1):
            es = lines[i].strip().split()
            dat_lst.append([])
            for v in range(self.n_vis):
                bin_val = int2binary(es[0][v])
                for u in range(2):
                    dat_lst[i-1].append( int(bin_val[u]) )
            self.prob_states[i-1] = float(es[1])
            #if self.prob_states[i-1] < 1e-10:
            #    self.prob_states[i-1] = 1e-10

        dat_lst = nd.array(dat_lst)

        #self.enum_states = nd.one_hot(dat_lst, self.n_val).reshape([-1, self.n_vis * self.n_val]).copyto(self.ctx)
        self.enum_states = dat_lst.copyto(self.ctx)
        sys.stderr.write("Exact states info loaded!\n")
        return

    def load_enum_RGs(self, fn):
        #use gsl-histogram result [0, 0.5)
        lines = open(fn, 'r').readlines()
        n_states = int(lines[0].split()[0])

        self.prob_RGs = nd.zeros([n_states], ctx=self.ctx)

        for i in range(1, n_states+1):
            es = lines[i].strip().split()
            self.prob_RGs[i-1] = float(es[1]) + 1e-10 #log(1/20000)
            assert( int(float(es[0])*2) == i-1 )

        sys.stderr.write("Rg info from obs data loaded!\n")
        return

    def check_KL(self):
        #print(self.enum_states)
        ph_act = nd.dot(self.enum_states, self.W) + self.hb
        vt = nd.dot(self.enum_states, self.vb)
        ht = nd.sum(-nd.log(nd.sigmoid(-ph_act)), axis=1)
        p_th = nd.softmax(vt+ht) + 1e-10
        KL = nd.sum(self.prob_states * nd.log(self.prob_states/p_th+1e-10))
        return KL.asnumpy()[0]

    def check_rgKL(self, v):
        tmp = v.reshape([-1, 2])
        ndx = tmp[:,0] * 2 + tmp[:,1]
        pv = np.ones(shape=self.prob_RGs.shape)
        for intcoor in ndx.asnumpy().reshape([-1, self.n_vis]):
            rg2 = cal_rg2(int2xy(intcoor))
            i = int(rg2*2)
            if i>=len(pv): i=len(pv)-1
            pv[i] += 1
        pv /= np.sum(pv)
        prg = self.prob_RGs.asnumpy()
        KL = np.sum(prg * np.log(prg/pv))
        return KL

    def sample_h_given_v(self, v0):
        h1_prob = self.propup(v0)
        h1 = h1_prob > nd.random_uniform(low=0.0, high=1.0, shape=h1_prob.shape, ctx=self.ctx)
        return [h1_prob, h1]

    def sample_v_given_h_without_softmax(self, h0):
        v1_prob = self.propdown(h0)
        #!!! what !!!
        v1_prob = nd.sigmoid(v1_prob)
        v1 = v1_prob > nd.random_uniform(low=0.0, high=1.0, shape=v1_prob.shape, ctx=self.ctx)
        return [v1_prob, v1]

    def propup(self, v):
        pre_sigmoid_activation = nd.dot(v, self.W) + self.hb
        return nd.sigmoid(pre_sigmoid_activation)

    def propdown(self, h):
        return nd.dot(h, self.W.T) + self.vb

    def gibbs_hvh(self, h0):
        v1_prob, v1 = self.sample_v_given_h_without_softmax(h0)
        h1_prob, h1 = self.sample_h_given_v(v1)
        return [v1_prob, v1, h1_prob, h1]

    def get_free_energy(self, v):
        x = nd.dot( v, self.W ) + self.hb
        vt = nd.dot( v, self.vb )
        ht = nd.sum( nd.log( 1.0 + nd.exp(x) ), axis=1 )
        fe = -ht-vt #free energy, how to prevent scale
        return nd.mean( fe )

    def reconstruct(self, v):
        h = nd.sigmoid( nd.dot( v, self.W ) + self.hb )
        reconstructed_v_prob = nd.sigmoid( nd.dot( h, self.W.T ) + self.vb )
        return reconstructed_v_prob

    def yeardream(self, v, n, k=1):
        sys.stderr.write("Year dreaming ...\n")
        n_sample = v.shape[0]
        #for i in range(20):
        #    h_prob, h0 = self.sample_h_given_v(v)
        #    v_prob, v = self.sample_v_given_h_without_softmax(h0)
        for i in range(n):
            for cdk in range(k):
                h_prob, h0 = self.sample_h_given_v(v)
                v_prob, v = self.sample_v_given_h_without_softmax(h0)
            #ndx = nd.argmax(v.reshape([-1,self.n_val]), axis=1)
            tmp = v.reshape([-1, 2])
            ndx = tmp[:,0] * 2 + tmp[:,1]
            for j, k in enumerate(ndx.asnumpy()):
                sys.stdout.write(str(int(k)))
                if j % self.n_vis == (self.n_vis-1):
                    sys.stdout.write("\n")
        return

    def load(self, fn):
        lines = open(fn, 'r').readlines()
        elems = lines[-1].split()
        last_epoch = int(elems[0])
        pos = 1
        for i in range(self.n_node):
            for j in range(self.n_hid):
                self.W[i,j] = float(elems[pos])
                pos += 1
        for i in range(self.n_node):
            self.vb[i] = float(elems[pos])
            pos += 1
        for j in range(self.n_hid):
            self.hb[j] = float(elems[pos])
            pos += 1

        sys.stderr.write("Loading weights and restart from epoch=%d\n" % last_epoch)
        return last_epoch

    def save(self, fn, epoch=0):
        with open(fn, 'a') as fp:
            fp.write("%d " % epoch)
            for i in range(self.n_node):
                for j in range(self.n_hid):
                    fp.write("%6.4f " % self.W.asnumpy()[i,j])
            for i in range(self.n_node):
                fp.write("%6.4f " % self.vb.asnumpy()[i])
            for j in range(self.n_hid):
                fp.write("%6.4f " % self.hb.asnumpy()[j])
            fp.write("\n")
        return

def get_data( fn, n_vis, n_val ):
    if fn.isdigit():
        #random
        num_data = int(fn)
        prob = nd.ones([num_data * n_vis, n_val]) / n_val
        dat_lst = nd.sample_multinomial(prob)
        sys.stderr.write("Generating random data: nv= %d, nd= %d\n" % (n_vis, num_data))
    else:
        #read from file
        with open(fn, 'r') as fp:
            lines = fp.readlines()
            es = lines[0].split()
            nl = int(es[0])
            nv = int(es[1])
            dat_lst = []
            sys.stderr.write("Loading data: nv= %d, nd= %d\n" % (nv, nl))
            for l in range(1, nl+1):
                dat_lst.append([])
                for i in range(nv):
                    bin_val = int2binary(lines[l][i])
                    for j in range(2):
                        dat_lst[l-1].append(int(bin_val[j]))
            data = nd.array(dat_lst)
    #data = nd.one_hot(dat_lst, n_val).reshape([-1, n_vis * n_val])
    return data

def run_rbm():
    #setup options
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", '--learning_rate', type=float, default=0.1, help="learning rate")
    parser.add_argument("-k", '--cdk', type=int, default=1, help="CD-k steps")
    parser.add_argument("-d", '--data', type=str, default="0", help="number/file of data")
    parser.add_argument("-t", '--train', type=int, default=0, help="training steps")
    parser.add_argument("-g", '--generate', type=int, default=0, help="generate steps")
    parser.add_argument("-b", '--batch_size', type=int, default=0, help="batch size")
    parser.add_argument("-v", '--visible', type=int, default=1, help="visible node")
    parser.add_argument("-n", '--hidden', type=int, default=1, help="hidden node")
    parser.add_argument("-l", '--level', type=int, default=2, help="level of vis node")
    parser.add_argument("-w", '--weight', type=str, help="R/W weights parameter")
    parser.add_argument("-i", '--interval', type=int, default=10, help="interval for output")
    parser.add_argument("-s", '--rescaleL', type=float, default=0.0, help="L1(<0)/L2(>0) type weights decay")
    parser.add_argument("-m", '--momentum', type=float, default=0.0, help="momentum coefficient")
    parser.add_argument("-x", '--gpu', type=int, default=-1, help="using GPU x or CPU(-1)")
    parser.add_argument("-kl", '--KL', type=str, help="check KL")
    parser.add_argument("-rg", '--rgKL', type=str, help="check KL of RG")

    args = parser.parse_args()
    learning_rate = args.learning_rate
    training_epochs = args.train
    gen_steps = args.generate
    bs = args.batch_size
    k = args.cdk
    num_vis = args.visible
    num_hid = args.hidden
    num_val = args.level

    if args.gpu<0:
        ctx = mx.cpu()
        sys.stderr.write("Using CPUs\n")
    else:
        ctx = mx.gpu(args.gpu)
        sys.stderr.write("Using GPU %d\n" % args.gpu)

    if args.data != "0":
        data = get_data(args.data, num_vis, num_val)
        if args.gpu>=0: data = data.copyto(ctx)
    else:
        data = None

    # construct RBM
    sys.stderr.write("Creating RBM: nv= %d, nl= %d, nh= %d\n" % (num_vis, num_val, num_hid))
    rbm = RBM(n_vis=num_vis, n_hid=num_hid, n_val=num_val, ctx=ctx)
    rbm.L_coeff = args.rescaleL
    rbm.M_coeff = args.momentum

    # load weights
    epoch_start = 0
    if args.weight is not None:
        if os.path.isfile(args.weight):
            epoch_start = rbm.load(args.weight)

    # load KL reference
    if args.KL is not None: rbm.load_enum_states(args.KL)
    if args.rgKL is not None: rbm.load_enum_RGs(args.rgKL)

    # train
    if training_epochs > 0 and data is not None:
        if epoch_start == 0: rbm.check_status(data, 0)
        for epoch in range(epoch_start, training_epochs):
            rbm.contrastive_divergence(input=data, lr=learning_rate, cdk=k,
                batch_size=bs, shuffle=(epoch>0))
            if (epoch+1) % args.interval == 0:
                rbm.check_status(data, epoch+1)
                sys.stdout.flush()
                sys.stderr.flush()
                #save weights (append), in one line
                if args.weight is not None: rbm.save(args.weight, epoch+1)

    if gen_steps > 0:
        if data is not None:
            rbm.yeardream(data, gen_steps, k)

if __name__ == "__main__":
    run_rbm()
