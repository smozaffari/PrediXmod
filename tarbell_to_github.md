###Directions for Accessing GitHub from a Secure Server 

I figured out how to _push_ to GitHub from a **secure** cluster like cri tarbell.

Although we can _pull_, _commit_, _branch_, _remote sync_, and _clone_ (i.e. do mostly everything) from a workhorse node 
like **lmem**, we cannot _push_ to GitHub from there. We indeed need to use a configured node like **cri-syncmon.cri.uchicago.edu**. 
However, since this is through ssh, the following steps are also necessary:

1. You need to create an [ssh key and passphrase](link1) 
2. You need to add the key to your GitHub account (see link above)
3. When cloning a GitHub repository, be sure to [clone in ssh](link2) 

You can now **_push_** your local commits on **tarbell** to your GitHub repository. I will be primarily working on (and committing from) 
the lmem node but will use the cri-syncmon node for pushing changes to GitHub at the very end of the workday. 
I'm happy to help anyone with questions on this process. 

[link1]: https://help.github.com/articles/generating-ssh-keys/
[link2]: https://help.github.com/articles/which-remote-url-should-i-use/

