<!--
{
  "webtitle": "Modules \u2014 specbvp documentation",
  "codeblocks": false
}
-->

# specbvp.polybases

Polynomial bases.

## PolyBasis

<pre class="py-sign"><b><em>class</em></b> specbvp.polybases.<b>PolyBasis</b>()</pre>

ABC for a polynomial basis.

- Returns *operators* derived from `PolyOpr`.
- Returns *node sets* derived from `NodeSet`.

### polys()

<pre class="py-sign">PolyBasis.<b>polys</b>(<em>self</em>)</pre>

Returns an *operator* that evaluates polynomials.

### derivs()

<pre class="py-sign">PolyBasis.<b>derivs</b>(<em>self</em>, order)</pre>

Returns an *operator* that differentiates polynomials.

<b>Parameters</b>

<p><span class="vardef"><code>order</code> : <em>int</em></span></p>

<dl><dd>
  Order of the desired derivative starting from one.
</dd></dl>

### integax()

<pre class="py-sign">PolyBasis.<b>integax</b>(<em>self</em>, weighted=<span>False</span>)</pre>

Returns an *operator* that integrates polynomials over `[a,x]`.

<b>Parameters</b>

<p><span class="vardef"><code>weighted</code> : <em>bool = False</em></span></p>

<dl><dd>
  Polynomials are multiplied by <code>x</code>, if <em>True</em>.
</dd></dl>

### integxb()

<pre class="py-sign">PolyBasis.<b>integxb</b>(<em>self</em>, weighted=<span>False</span>)</pre>

Returns an *operator* that integrates polynomials over `[x,b]`.

<b>Parameters</b>

<p><span class="vardef"><code>weighted</code> : <em>bool = False</em></span></p>

<dl><dd>
  Polynomials are multiplied by <code>x</code>, if <em>True</em>.
</dd></dl>

### nodes()

<pre class="py-sign">PolyBasis.<b>nodes</b>(<em>self</em>) → <em>dict</em></pre>

Returns a dictionary with the available node sets.

## PolyOpr

<pre class="py-sign"><b><em>class</em></b> specbvp.polybases.<b>PolyOpr</b>()</pre>

ABC for operators on a polynomial sequence.

### setnodes()

<pre class="py-sign">PolyOpr.<b>setnodes</b>(<em>self</em>, nodes)</pre>

Defines the output points and returns the instance.

<b>Parameters</b>

<p><span class="vardef"><code>nodes</code> : <em>number | array-like</em></span></p>

<dl><dd>
  Collocation point(s) within <code>[a,b]</code>.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>self</em></span></p>

<dl><dd>
  The instance itself.
</dd></dl>

### setpolys()

<pre class="py-sign">PolyOpr.<b>setpolys</b>(<em>self</em>, *indices)</pre>

Defines a polynomial sequence and returns the instance.

<b>Parameters</b>

<p><span class="vardef"><code>indices</code> : <em>*int</em></span></p>

<dl><dd>
  Indices of the polynomials to include.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>self</em></span></p>

<dl><dd>
  The instance itself.
</dd></dl>

### asdict()

<pre class="py-sign">PolyOpr.<b>asdict</b>(<em>self</em>) → <em>dict</em></pre>

Realizes the operator as an index-to-output mapping.

<b>Returns</b>

<p><span class="vardef"><em>dict</em></span></p>

<dl><dd>
  Maps the indices of polynomials to the output values.
</dd></dl>

### asmat()

<pre class="py-sign">PolyOpr.<b>asmat</b>(<em>self</em>)</pre>

Realizes the operator as a Vandermonde-like matrix.

<b>Returns</b>

<p><span class="vardef"><em>ndarray</em></span></p>

<dl><dd>
  The operator as a Vandermonde-like matrix (a).
</dd></dl>

(a) Columns are images of the polynomials tabulated at the nodes.

## NodeSet

<pre class="py-sign"><b><em>class</em></b> specbvp.polybases.<b>NodeSet</b>()</pre>

Set of nodes associated with a polynomials basis.

- Represents a set of distinct points in `[a,b]`.
- May hold weights of the accompanying quadrature rule. 

<b>Attributes</b>

<p><span class="vardef"><code>nodes</code> : <em>ndarray</em></span></p>

<dl><dd>
  Points as a flat numpy array.
</dd></dl>

<p><span class="vardef"><code>weights</code> : <em>ndarray | None</em></span></p>

<dl><dd>
  Weights as a flat numpy array, if any.
</dd></dl>

### setnum()

<pre class="py-sign">NodeSet.<b>setnum</b>(<em>self</em>, number)</pre>

Defines the number of points and returns the instance.

<b>Parameters</b>

<p><span class="vardef"><code>size</code> : <em>int</em></span></p>

<dl><dd>
  Number of nodes.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>self</em></span></p>

<dl><dd>
  The instance itself.
</dd></dl>

## Legendre

<pre class="py-sign"><b><em>class</em></b> specbvp.polybases.<b>Legendre</b>()</pre>

Basis formed by the Legendre polynomials.

## Chebyshev

<pre class="py-sign"><b><em>class</em></b> specbvp.polybases.<b>Chebyshev</b>()</pre>

Basis formed by the Chebyshev polynomials of the 1st kind.