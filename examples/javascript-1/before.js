const person = {
    "name": "Alice",
    "age": 28,
    "email": "alice@example.com",
    "hobbies": ["reading", "cycling", "coding"]
};

console.log("Name:", person.name);
console.log("Age:", person.age);
console.log("Hobbies:", person.hobbies.join(", "));

person.age += 1;
console.log("Updated Age:", person.age);

const jsonString = JSON.stringify(person);
console.log("JSON String:", jsonString);
