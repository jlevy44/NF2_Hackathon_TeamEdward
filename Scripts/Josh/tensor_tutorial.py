import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'
import tensorflow as tf

node1 = tf.constant(3.0, dtype=tf.float32)
node2 = tf.constant(4.0) # also tf.float32 implicitly
print(node1,node2)
sess = tf.Session()
print(sess.run([node1, node2]))

node3 = tf.add(node1,node2)
print("node3: ", node3)
print("sess.run(node3): ",sess.run(node3))

a = tf.placeholder(tf.float32)
b = tf.placeholder(tf.float32)
adder_node = a + b  # + provides a shortcut for tf.add(a, b)

print(sess.run(adder_node, {a: 3, b:4.5}))
print(sess.run(adder_node, {a: [1,3], b: [2, 4]}))
