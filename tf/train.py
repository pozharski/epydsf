from argparse import ArgumentParser, RawDescriptionHelpFormatter

def main():
    headerhelp = \
'''
Train NN to recognize melting points.

'''
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
    parser.add_argument('--t1', type=float, default=10.0, help='Lower temperature limit, degrees Celsius.')
    parser.add_argument('--t2', type=float, default=90.0, help='Upper temperature limit, degrees Celsius.')
    parser.add_argument('--ht', type=float, default=1.0, help='Temperature step.')
    parser.add_argument('--tmrange', default='30,70,1', help='Comma-separated list defining training Tm range')
    parser.add_argument('--dtrange', default='1,10,0.5', help='Comma-separated list defining training dT range')

    args = parser.parse_args()

    from data import baseline, melt
    from scipy import arange, array
    import random
    random.seed()
    tm1,tm2,dtm = [float(x) for x in args.tmrange.split(',')]
    dt1,dt2,ddt = [float(x) for x in args.dtrange.split(',')]
    
    t = arange(args.t1,args.t2,args.ht)
    samples = array([[melt([tm,dt]).values(t).tolist()+[tm] for dt in arange(dt1,dt2,ddt)] for tm in arange(tm1,tm2,dtm)])
    samples = samples.reshape((samples.shape[0]*samples.shape[1],samples.shape[2]))
    print("Training set contains %d curves" % samples.shape[0])
    tms = samples[:,-1]
    samples = samples[:,:-1]
    
    input_num_units = samples.shape[1]
    hidden_num_units = 512
    output_num_units = samples.shape[1]
    
    x = tf.placeholder(tf.float32, [None, input_num_units])
    y = tf.placeholder(tf.float32, [None, output_num_units])

if __name__ == "__main__":
    main()
