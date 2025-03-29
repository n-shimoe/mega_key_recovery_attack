# Proof of Concept (PoC) for RSA key recovery attacks on Unpatched MEGA
This code verifies our novel RSA key recovery attack on unpatched MEGA with 4 queries in the simulated and simplified setting.
This code also includes [Heninger-Ryan's RSA key recovery attack with 6 queries](https://eprint.iacr.org/2022/914) and the variant of our attack using [Lu-Zhang-Peng-Lin's method](https://link.springer.com/chapter/10.1007/978-3-662-48797-6_9).

## Requirements
SageMath with Python 3.11.1. SageMath 9.8 is recommended.

Note that our implementation may produce the unintended result in the older version of Python or SageMath. 

## Code Organization
<dl>
    <dt>poc_optimized.sage</dt>
    <dd>our novel RSA key recovery attack with 4 queries</dd>
    <dt>poc_heninger_ryan.sage</dt>
    <dd>Heninger-Ryan's RSA key recovery attack with 6 queries</dd>
    <dt>poc_lzpl_variant.sage</dt>
    <dd>the variant of our attack using Lu-Zhang-Peng-Lin's method</dd>
</dl>