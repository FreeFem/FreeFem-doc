See [Loop example](../examples/#loop).

## for
For loop.
```freefem
for (int i = 0; i < N; ++i){
	...
}
```

## if
If condition.
```freefem
if (condition){
	...
}
else{
	...
}
```

## else
See [if](#if).

## while
While loop.
```freefem
while (condition){
	...
}
```

## continue
Continue a loop.
```freefem
for (int i = 0; i < N; ++i){
	...
	if (condition) continue;
	...
}
```

## break
Break a loop.
```freefem
while (condition1){
	...
	if (condition) break;
	...
}
```

## try
Try a part of code.
```freefem
try{
	...
}
catch(...){
	...
}
```

See [Basic error handling example](../examples/#basic-error-handling) and [Error handling example](../examples/#error-handling).

## catch
Catch an error, see [try](#try)

## Implicit loop

Array with one index:
```freefem
for [i, ai : a]
```
If `:::freefem real[int] a(10)`, then `i=0:9` and `ai` is a reference to `a[i]`.

Array with two indices or matrix:
```freefem
for [i, j, aij : a]
```
If `:::freefem real[int] a(10, 11)`, then `i=0:9`, `j=1:10` and `aij` is a reference to `a(i, j)`.

See  [Implicit loop example](../examples/implicit-loop).
