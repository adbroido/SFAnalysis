import urllib2
from multiprocessing.dummy import Pool as ThreadPool
import time


"""
 based on http://chriskiehl.com/article/parallelism-in-one-line/
"""

urls = [
  'http://www.python.org',
  'http://www.python.org/about/',
  'http://www.onlamp.com/pub/a/python/2003/04/17/metaclasses.html',
  'http://www.python.org/doc/',
  'http://www.python.org/download/',
  'http://www.python.org/getit/',
  'http://www.python.org/community/',
  'https://wiki.python.org/moin/',
  'http://planet.python.org/',
  'https://wiki.python.org/moin/LocalUserGroups',
  'http://www.python.org/psf/',
  'http://docs.python.org/devguide/',
  'http://www.python.org/community/awards/'
  # etc..
  ]

startfor = time.time()
results = []
for url in urls:
    result = urllib2.urlopen(url)
    results.append(result)
endfor = time.time()

# # ------- 4 Pool ------- #
# Make the Pool of workers
start4 = time.time()
pool = ThreadPool(4)
# Open the urls in their own threads
# and return the results
results = pool.map(urllib2.urlopen, urls)
#close the pool and wait for the work to finish
pool.close()
pool.join()
end4 = time.time()

# # ------- 8 Pool ------- #
# Make the Pool of workers
start8 = time.time()
pool = ThreadPool(8)
# Open the urls in their own threads
# and return the results
results = pool.map(urllib2.urlopen, urls)
#close the pool and wait for the work to finish
pool.close()
pool.join()
end8 = time.time()


# # ------- 13 Pool ------- #
# Make the Pool of workers
start13 = time.time()
pool = ThreadPool(13)
# Open the urls in their own threads
# and return the results
results = pool.map(urllib2.urlopen, urls)
#close the pool and wait for the work to finish
pool.close()
pool.join()
end13 = time.time()


print 'Single thread: %0.2f seconds' %(endfor-startfor)
print '4 Pool: %0.2f seconds' %(end4-start4)
print '8 Pool: %0.2f seconds' %(end8-start8)
print '13 Pool: %0.2f seconds' %(end13-start13)
