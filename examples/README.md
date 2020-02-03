To fetch and build XCFun as a dependency
----------------------------------------

```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test
```


To use an existing XCFun installation
-------------------------------------

```
$ mkdir build
$ cd build
$ cmake .. -DXCFun_DIR=<path/to/xcfun>/share/cmake/XCFun
$ make
$ make test
```
